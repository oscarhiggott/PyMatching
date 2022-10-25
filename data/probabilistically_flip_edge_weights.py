import re
import random
import typer


def invert_dem_probabilities_with_flip_probability(dem_fn: str, prob_invert: float) -> None:
    with open(dem_fn, "r") as f:
        lines = f.readlines()

    for line in lines:
        m = re.match(".*?error\((.+)\).+", line)
        if m is not None:
            p_str = m.group(1)
            p = float(p_str)
            if random.random() < prob_invert:
                line = re.sub(p_str, str(1 - p), line)
        print(line, end="")


if __name__ == "__main__":
    typer.run(invert_dem_probabilities_with_flip_probability)
