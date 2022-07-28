#include <gtest/gtest.h>
#include "alternating_tree.h"
#include "graph_fill_region.h"


std::vector<pm::GraphFillRegion>& getGfrs(){
    static std::vector<pm::GraphFillRegion> gfrs(15);
    return gfrs;
}

pm::Graph& getGraph() {
    static pm::Graph g(15);
    return g;
}


//pm::AltTreeEdge t(
//        const std::vector<pm::AltTreeEdge>& children,
//        size_t inner_region_id,
//        size_t outer_region_id,
//        bool root,
//        size_t io_loc_from_id,
//        size_t io_loc_to_id,
//        size_t parent_loc_from_id,
//        size_t parent_loc_to_id
//        ) {
//    pm::AltTreeNode* node;
//    pm::CompressedEdge parent_ce;
//    if (root) {
//        node = new pm::AltTreeNode(&getGfrs()[outer_region_id]);
//    } else {
//        node = new pm::AltTreeNode(
//                &getGfrs()[inner_region_id],
//                &getGfrs()[outer_region_id],
//                {
//                    &getGraph().nodes[io_loc_from_id],
//                    &getGraph().nodes[io_loc_to_id],
//                    0
//                }
//                );
//        parent_ce = {
//                &getGraph().nodes[parent_loc_from_id],
//                &getGraph().nodes[parent_loc_to_id],
//                0
//        };
//    }
//    auto edge = pm::AltTreeEdge(node, parent_ce);
//    for (auto& child : children) {
//        edge.alt_tree_node->add_child(child);
//    }
//    return edge;
//}

pm::AltTreeEdge t(
        std::vector<pm::AltTreeEdge> children,
        size_t inner_region_id,
        size_t outer_region_id,
        bool root = false
) {
    pm::AltTreeNode* node;
    pm::CompressedEdge parent_ce;
    if (root) {
        node = new pm::AltTreeNode(&getGfrs()[outer_region_id]);
    } else {
        node = new pm::AltTreeNode(
                &getGfrs()[inner_region_id],
                &getGfrs()[outer_region_id],
                {
                        &getGraph().nodes[inner_region_id],
                        &getGraph().nodes[outer_region_id],
                        0
                }
        );
    }
    auto edge = pm::AltTreeEdge(node, parent_ce);
    for (auto& child : children) {
        child.edge.loc_from = &getGraph().nodes[outer_region_id];
        child.edge.loc_to = child.alt_tree_node->inner_to_outer_edge.loc_from;
        edge.alt_tree_node->add_child(child);
    }
    return edge;
}


pm::AltTreeEdge example_tree() {
    return t(
            {
                    t({}, 7, 8
                    ),
                    t(
                            {
                                    t({}, 3, 4),
                                    t({}, 5, 6)
                            },
                            1, 2
                    )
            },
            -1, 0, true
    );
}

pm::AltTreeEdge example_tree_four_children() {
    return t(
            {
                    t({}, 1, 2),
                    t({}, 3, 4),
                    t({}, 5, 6),
                    t({}, 7, 8)
            },
            -1, 0, true);
}


void delete_alternating_tree(pm::AltTreeNode* root) {
    for (auto child: root->children) {
        delete_alternating_tree(child.alt_tree_node);
    }
    delete root;
}


TEST(AlternatingTree, FindRoot) {
    auto n = example_tree();
    ASSERT_EQ(
            n.alt_tree_node->children[1].alt_tree_node->children[0].alt_tree_node->find_root(), n.alt_tree_node);
    ASSERT_EQ(n.alt_tree_node->find_root(), n.alt_tree_node);
    delete_alternating_tree(n.alt_tree_node);
}


TEST(AlternatingTree, UnstableEraseInt) {
    std::vector<int> v = {6, 2, 4, 7, 9, 10};
    bool b1 = pm::unstable_erase(v, [](int x){return x==7;});
    ASSERT_EQ(v, std::vector<int>({6, 2, 4, 10, 9}));
    ASSERT_TRUE(b1);
    bool b2 = pm::unstable_erase(v, [](int x){return x == 6;});
    ASSERT_EQ(v, std::vector<int>({9, 2, 4, 10}));
    ASSERT_TRUE(b2);
    std::vector<int> w = {8};
    bool b3 = pm::unstable_erase(w, [](int x){return x == 8;});
    ASSERT_EQ(w, std::vector<int>({}));
    ASSERT_TRUE(b3);
    bool b4 = pm::unstable_erase(w, [](int x){return x == 0;});
    ASSERT_EQ(w, std::vector<int>({}));
    ASSERT_FALSE(b4);
}


TEST(AlternatingTree, UnstableEraseAltTreeEdge) {
    pm::AltTreeEdge x = example_tree_four_children();
    std::vector<pm::AltTreeEdge> xc = x.alt_tree_node->children;
    auto xc_copy = xc;
    ASSERT_EQ(xc, x.alt_tree_node->children);
    pm::unstable_erase(xc, [&xc_copy](pm::AltTreeEdge y){
        return y.alt_tree_node == xc_copy[1].alt_tree_node;});
    ASSERT_EQ(xc, std::vector<pm::AltTreeEdge>({xc_copy[0], xc_copy[3], xc_copy[2]}));
    pm::unstable_erase(xc, [&xc_copy](pm::AltTreeEdge y){
        return y.alt_tree_node == xc_copy[0].alt_tree_node;});
    ASSERT_EQ(xc, std::vector<pm::AltTreeEdge>({xc_copy[2], xc_copy[3]}));
    delete_alternating_tree(x.alt_tree_node);
}

TEST(AlternatingTree, AllNodesInTree) {
    pm::AltTreeEdge x = example_tree();
    ASSERT_EQ(x.alt_tree_node->all_nodes_in_tree(),
              std::vector<pm::AltTreeNode*>({
                  x.alt_tree_node,
                  x.alt_tree_node->children[1].alt_tree_node,
                  x.alt_tree_node->children[1].alt_tree_node->children[1].alt_tree_node,
                  x.alt_tree_node->children[1].alt_tree_node->children[0].alt_tree_node,
                  x.alt_tree_node->children[0].alt_tree_node
              }));
    delete_alternating_tree(x.alt_tree_node);
}

TEST(AlternatingTree, TreeEqual) {
    pm::AltTreeEdge x = example_tree();
    pm::AltTreeEdge x2 = example_tree();
    pm::AltTreeEdge y = example_tree_four_children();
    pm::AltTreeEdge y2 = example_tree_four_children();
    ASSERT_TRUE(x.alt_tree_node->tree_equal(*x2.alt_tree_node));
    ASSERT_TRUE(*x.alt_tree_node == *x2.alt_tree_node);
    ASSERT_TRUE(y.alt_tree_node->tree_equal(*y.alt_tree_node));
    ASSERT_TRUE(y.alt_tree_node->tree_equal(*y2.alt_tree_node));
    ASSERT_FALSE(x.alt_tree_node->tree_equal(*y.alt_tree_node));
    ASSERT_FALSE(*x.alt_tree_node == *y.alt_tree_node);
    y2.alt_tree_node->children.pop_back();
    ASSERT_FALSE(y.alt_tree_node->tree_equal(*y2.alt_tree_node));
    x2.alt_tree_node->children[0].alt_tree_node->outer_region = nullptr;
    ASSERT_FALSE(*x.alt_tree_node == *x2.alt_tree_node);
    delete_alternating_tree(x.alt_tree_node);
    delete_alternating_tree(x2.alt_tree_node);
    delete_alternating_tree(y.alt_tree_node);
    delete_alternating_tree(y2.alt_tree_node);
}


TEST(AlternatingTree, BecomeRoot) {
    auto y1 = t(
            {t({}, 1, 2)},
            -1,
            0,
            true
            );
    auto y2 = t(
            {t({}, 1, 0)}, -1, 2, true
            );
    auto v = y1.alt_tree_node->children[0].alt_tree_node;
    v->become_root();
    ASSERT_EQ(*v, *y2.alt_tree_node);

    auto x1 = t(
            {
                t(
                        {
                            t({}, 10, 11),
                            t({}, 8, 9),
                            t({ t({}, 6, 7)}, 12, 13)
                        }, 4, 5),
                t({}, 2, 3)
            },
            -1,
            1,
            true);

    auto x2 = t(
            {
                t({}, 6, 7),
                t({
                            t({}, 10, 11),
                            t({}, 8, 9),
                            t({t({}, 2, 3)}, 4, 1)},
                12, 5
            )
            },
            -1, 13, true);
    auto c = x1.alt_tree_node->children[0].alt_tree_node->children[2].alt_tree_node;
    c->become_root();
    ASSERT_EQ(*c, *x2.alt_tree_node);
}
