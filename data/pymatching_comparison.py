import math
from typing import List, Callable, Set, Tuple

import pymatching
import networkx as nx
import numpy as np
import stim
import os
import typer
import tqdm
import time


def nx_graph_to_pymatching(g: nx.Graph) -> pymatching.Matching:
    for _, _, d in g.edges(data=True):
        d["fault_ids"] = set_bits(d["observables"])
    max_node_id = max(g.nodes())
    max_fault_id = max(max(d["fault_ids"], default=0) for _, _, d in g.edges(data=True))
    g.add_edge(max_node_id, max_node_id + 1, weight=1, fault_ids=max_fault_id)
    for i in range(len(g.nodes)):
        g.nodes[i]["is_boundary"] = g.nodes[i].get("boundary", False)
    g.nodes[max_node_id + 1]["is_boundary"] = True
    return pymatching.Matching(g)


def detector_error_model_to_discretised_pymatching_graph(
        model: stim.DetectorErrorModel, num_buckets: int = 1000
) -> pymatching.Matching:
    graph = detector_error_model_to_nx_graph(model)
    graph = discretize_weights(graph, num_buckets)
    return nx_graph_to_pymatching(graph)


def detector_error_model_to_pymatching_graph(
        model: stim.DetectorErrorModel
) -> pymatching.Matching:
    graph = detector_error_model_to_nx_graph(model)
    return nx_graph_to_pymatching(graph)


def int_to_binary_array(n: int, num_bits: int) -> np.ndarray:
    obs_list = [int(i) for i in bin(n)[:1:-1]]
    obs_list += [0] * (num_bits - len(obs_list))
    return np.array(obs_list, dtype=np.uint8)


def set_bits(n: int) -> Set[int]:
    return {i for i, c in enumerate(bin(n)[:1:-1]) if c == '1'}


def discretize_weights(g: nx.Graph, num_buckets) -> nx.Graph:
    max_weight = max(abs(float(e[2].get("weight", 1))) for e in g.edges(data=True))
    for u, v, d in g.edges(data=True):
        d["weight"] = 2 * round(num_buckets * (float(d.get("weight", 1)) / max_weight))
    return g


def detector_error_model_to_nx_graph(model: stim.DetectorErrorModel) -> nx.Graph:
    """Convert a stim error model into a NetworkX graph.
    From: https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb"""

    g = nx.Graph()
    boundary_node = model.num_detectors
    g.add_node(boundary_node, boundary=True, coords=[-1, -1, -1])

    def handle_error(p: float, dets: List[int], frame_changes: int):
        if p == 0:
            return
        if len(dets) == 0:
            # No symptoms for this error.
            # Code probably has distance 1.
            # Accept it and keep going, though of course decoding will probably perform terribly.
            return
        if len(dets) == 1:
            dets = [dets[0], boundary_node]
        if len(dets) > 2:
            raise NotImplementedError(
                f"Error with more than 2 symptoms can't become an edge or boundary edge: {dets!r}.")
        if g.has_edge(*dets):
            edge_data = g.get_edge_data(*dets)
            old_p = edge_data["error_probability"]
            old_frame_changes = edge_data["observables"]
            # If frame changes differ, the code has distance 2; just keep whichever was first.
            if old_frame_changes == frame_changes:
                p = p * (1 - old_p) + old_p * (1 - p)
                g.remove_edge(*dets)
        g.add_edge(*dets, weight=math.log((1 - p) / p), observables=frame_changes, error_probability=p)

    def handle_detector_coords(detector: int, coords: np.ndarray):
        g.add_node(detector, coords=coords)

    eval_model(model, handle_error, handle_detector_coords)

    return g


def eval_model(
        model: stim.DetectorErrorModel,
        handle_error: Callable[[float, List[int], int], None],
        handle_detector_coords: Callable[[int, np.ndarray], None]):
    """Interprets the error model instructions, taking care of loops and shifts.
    Adapted from: https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb

    Makes callbacks as error mechanisms are declared, and also when detector
    coordinate data is declared.
    """
    det_offset = 0
    coords_offset = np.zeros(100, dtype=np.float64)

    def _helper(m: stim.DetectorErrorModel, reps: int):
        nonlocal det_offset
        nonlocal coords_offset
        for _ in range(reps):
            for instruction in m:
                if isinstance(instruction, stim.DemRepeatBlock):
                    _helper(instruction.body_copy(), instruction.repeat_count)
                elif isinstance(instruction, stim.DemInstruction):
                    if instruction.type == "error":
                        dets: List[int] = []
                        frames: int = 0
                        t: stim.DemTarget
                        p = instruction.args_copy()[0]
                        for t in instruction.targets_copy():
                            if t.is_relative_detector_id():
                                dets.append(t.val + det_offset)
                            elif t.is_logical_observable_id():
                                frames ^= 2**t.val
                            elif t.is_separator():
                                # Treat each component of a decomposed error as an independent error.
                                # (Ideally we could configure some sort of correlated analysis; oh well.)
                                handle_error(p, dets, frames)
                                frames = 0
                                dets = []
                        # Handle last component.
                        handle_error(p, dets, frames)
                    elif instruction.type == "shift_detectors":
                        det_offset += instruction.targets_copy()[0]
                        a = np.array(instruction.args_copy())
                        coords_offset[:len(a)] += a
                    elif instruction.type == "detector":
                        a = np.array(instruction.args_copy())
                        for t in instruction.targets_copy():
                            handle_detector_coords(t.val + det_offset, a + coords_offset[:len(a)])
                    elif instruction.type == "logical_observable":
                        pass
                    else:
                        raise NotImplementedError()
                else:
                    raise NotImplementedError()
    _helper(model, 1)


def obs_array_to_obs_mask(obs_array: np.ndarray) -> int:
    obs_mask: int = 0
    for i, b in enumerate(obs_array):
        obs_mask ^= 1 << i if b else 0
    return obs_mask


def decode_with_pymatching(model: stim.DetectorErrorModel, samples: np.ndarray,
                           num_buckets: int = 1000, num_neighbours: int = -1,
                           verbose: bool = False, show_time: bool = False, show_failures: bool = False,
                           show_weights: bool = False, max_shots: int = None, show_predictions: bool = False) -> None:
    num_neighbours = None if num_neighbours == -1 else num_neighbours
    if num_buckets is not None:
        matching = detector_error_model_to_discretised_pymatching_graph(model, num_buckets=num_buckets)
    else:
        matching = detector_error_model_to_pymatching_graph(model)

    dets = samples[:, 0:model.num_detectors]
    obs = samples[:, model.num_detectors:model.num_detectors+model.num_observables]
    num_failures = 0
    it = tqdm.tqdm(dets) if verbose else dets
    if num_neighbours is None:
        matching.matching_graph.compute_all_pairs_shortest_paths()
    i = 0
    t0 = time.time()
    for s in it:
        obs_pred, weight = matching.decode(s, return_weight=True, num_neighbours=num_neighbours)
        if show_weights:
            if num_buckets is not None:
                print(int(weight))
            else:
                print(weight)
        if show_failures:
            num_failures += not np.array_equal(obs[i, :].astype(np.uint8), obs_pred.astype(np.uint8))
        if show_predictions:
            obs_mask = obs_array_to_obs_mask(obs_pred)
            print(obs_mask)
        i += 1
        if max_shots is not None and i >= max_shots:
            break
    t1 = time.time()
    if show_time:
        print(f"{(t1-t0)*1e6/i}")
    if show_failures:
        print(f"{num_failures} {i}")


def decode_with_pymatching_from_files(model_fn: str, samples_fn: str, num_buckets: int = None,
                                      num_neighbours: int = -1, verbose: bool = False,
                                      show_time: bool = False, show_failures: bool = False,
                                      show_weights: bool = False, max_shots: int = None,
                                      show_predictions: bool = False):
    dem = stim.DetectorErrorModel.from_file(model_fn)
    _, extension = os.path.splitext(samples_fn)
    samples = stim.read_shot_data_file(path=samples_fn, format=extension[1:],
                                       num_measurements=0, num_detectors=dem.num_detectors,
                                       num_observables=dem.num_observables)
    decode_with_pymatching(dem, samples, num_buckets=num_buckets,
                           num_neighbours=num_neighbours, verbose=verbose,
                           show_time=show_time, show_failures=show_failures, show_weights=show_weights,
                           max_shots=max_shots, show_predictions=show_predictions)


if __name__ == "__main__":
    # This command line tool can be used to generate predictions and MWPM solution weights
    # using the previous version of pymatching (e.g. pymatching v0.7.0), for testing purposes.
    # To run this script, first `pip install pymatching==0.7.0 networkx numpy stim typer tqdm`
    typer.run(decode_with_pymatching_from_files)
