//
// Created by Nikolay Yakovets on 2018-02-02.
//

#include "SimpleEstimator.h"
#include "SimpleEvaluator.h"

using namespace std;

SimpleEvaluator::SimpleEvaluator(std::shared_ptr<SimpleGraph> &g) {

    // works only with SimpleGraph
    graph = g;
    est = nullptr; // estimator not attached by default
}

void SimpleEvaluator::attachEstimator(std::shared_ptr<SimpleEstimator> &e) {
    est = e;
}

void SimpleEvaluator::prepare() {
    // if attached, prepare the estimator
    if(est != nullptr) est->prepare();

    // prepare other things here.., if necessary

    //createAggregateIndex();
    createExhaustiveIndex();
}

void SimpleEvaluator::createExhaustiveIndex() {
    // exhaustive indexes: SOP, PSO, POS, OSP
    exh_indexes.POS.resize(graph->getNoLabels());
    exh_indexes.PSO.resize(graph->getNoLabels());
    for(uint32_t j = 0; j < graph->getNoVertices(); j++) {
        for (auto edge: graph->adj[j]) {
            //edge.first = edge type, edge.second = out node, j = in node
            // POS = edge type -> (out node, in node)
            exh_indexes.POS[edge.first].push_back(std::make_pair(edge.second, j));
            // PSO = edge type -> (in node, out node)
            exh_indexes.PSO[edge.first].push_back(std::make_pair(j, edge.second));
        }
    }
}

void SimpleEvaluator::createAggregateIndex() {
    // aggregate indexes: SP, SO, PS, PO, OS, OP, S, P, and O
    agg_indexes.SP.resize(graph->getNoVertices());
    agg_indexes.SO.resize(graph->getNoVertices());
    agg_indexes.PS.resize(graph->getNoLabels());
    agg_indexes.PO.resize(graph->getNoLabels());
    agg_indexes.OS.resize(graph->getNoVertices());
    agg_indexes.OP.resize(graph->getNoVertices());
    for(uint32_t j = 0; j < graph->getNoVertices(); j++) {
        for (auto edge: graph->adj[j]) {
            //edge.first = edge type, edge.second = out node, j = in node
            // SP = in node -> edge type
            agg_indexes.SP[j].push_back(edge.first);
            // SO = in node -> out node
            agg_indexes.SO[j].push_back(edge.second);
            // PS = edge type -> in node
            agg_indexes.PS[edge.first].push_back(j);
            // PO = edge type -> out node
            agg_indexes.PO[edge.first].push_back(edge.second);
            // OS = out node -> in node
            agg_indexes.OS[edge.second].push_back(j);
            // OP = out node + edge type
            agg_indexes.OP[edge.second].push_back(edge.first);
        }
    }
}

cardStat SimpleEvaluator::computeStats(std::shared_ptr<SimpleGraph> &g) {

    cardStat stats {};

    for(int source = 0; source < g->getNoVertices(); source++) {
        if(!g->adj[source].empty()) stats.noOut++;
    }

    stats.noPaths = g->getNoDistinctEdges();

    for(int target = 0; target < g->getNoVertices(); target++) {
        if(!g->reverse_adj[target].empty()) stats.noIn++;
    }

    return stats;
}

std::shared_ptr<SimpleGraph> SimpleEvaluator::project(uint32_t projectLabel, bool inverse, std::shared_ptr<SimpleGraph> &in) {

    auto out = std::make_shared<SimpleGraph>(in->getNoVertices());
    out->setNoLabels(in->getNoLabels());

    if(!inverse) {
        // going forward
        for(uint32_t source = 0; source < in->getNoVertices(); source++) {
            for (auto labelTarget : in->adj[source]) {

                auto label = labelTarget.first;
                auto target = labelTarget.second;

                if (label == projectLabel)
                    out->addEdge(source, target, label);
            }
        }
    } else {
        // going backward
        for(uint32_t source = 0; source < in->getNoVertices(); source++) {
            for (auto labelTarget : in->reverse_adj[source]) {

                auto label = labelTarget.first;
                auto target = labelTarget.second;

                if (label == projectLabel)
                    out->addEdge(source, target, label);
            }
        }
    }

    return out;
}

std::shared_ptr<SimpleGraph> SimpleEvaluator::project_exh_index(uint32_t projectLabel, bool inverse, std::shared_ptr<SimpleGraph> &in) {
    auto out = std::make_shared<SimpleGraph>(in->getNoVertices());
    out->setNoLabels(in->getNoLabels());

    auto PSO = exh_indexes.PSO[projectLabel];
    auto POS = exh_indexes.POS[projectLabel];

    if (!inverse) {
        // forward
        for (auto edge : PSO) {
            // edge.first = in node, edge.second = out node
            out->addEdge(edge.first, edge.second, projectLabel);
        }
    } else {
        // backward
        for (auto edge : POS) {
            // edge.first = out node, edge.second = in node
            out->addEdge(edge.first, edge.second, projectLabel);
        }
    }

    return out;
}

std::shared_ptr<SimpleGraph> SimpleEvaluator::project_agg_index(uint32_t projectLabel, bool inverse, std::shared_ptr<SimpleGraph> &in) {

    // find all nodes connected to edge type projectLabel
    std::vector<int> PS = agg_indexes.PS[projectLabel];
    std::vector<int> PO = agg_indexes.PO[projectLabel];

    auto out = std::make_shared<SimpleGraph>(in->getNoVertices());
    out->setNoLabels(in->getNoLabels());

    if (!inverse) {
        // forward
        for (auto inNode : PS) {
            for (auto outNode: PO) {
                out->addEdge(inNode, outNode, projectLabel);
            }
        }
    } else {
        // backward
        for (auto outNode : PO) {
            for (auto inNode: PS) {
                out->addEdge(inNode, outNode, projectLabel);
            }
        }
    }

    return out;
}

std::shared_ptr<SimpleGraph> SimpleEvaluator::join(std::shared_ptr<SimpleGraph> &left, std::shared_ptr<SimpleGraph> &right) {

    auto out = std::make_shared<SimpleGraph>(left->getNoVertices());
    out->setNoLabels(1);

    for(uint32_t leftSource = 0; leftSource < left->getNoVertices(); leftSource++) {
        for (auto labelTarget : left->adj[leftSource]) {

            int leftTarget = labelTarget.second;
            // try to join the left target with right source
            for (auto rightLabelTarget : right->adj[leftTarget]) {

                auto rightTarget = rightLabelTarget.second;
                out->addEdge(leftSource, rightTarget, 0);

            }
        }
    }

    return out;
}

//std::shared_ptr<SimpleGraph> SimpleEvaluator::evaluate_aux(RPQTree *q) {
//
//    // evaluate according to the AST bottom-up
//
//    if(q->isLeaf()) {
//        // project out the label in the AST
//        std::regex directLabel (R"((\d+)\+)");
//        std::regex inverseLabel (R"((\d+)\-)");
//
//        std::smatch matches;
//
//        uint32_t label;
//        bool inverse;
//
//        if(std::regex_search(q->data, matches, directLabel)) {
//            label = (uint32_t) std::stoul(matches[1]);
//            inverse = false;
//        } else if(std::regex_search(q->data, matches, inverseLabel)) {
//            label = (uint32_t) std::stoul(matches[1]);
//            inverse = true;
//        } else {
//            std::cerr << "Label parsing failed!" << std::endl;
//            return nullptr;
//        }
//
//        //return SimpleEvaluator::project(label, inverse, graph);
//        //return SimpleEvaluator::project_agg_index(label, inverse, graph);
//        return SimpleEvaluator::project_exh_index(label, inverse, graph);
//    }
//
//    if(q->isConcat()) {
//
//        // evaluate the children
//        auto leftGraph = SimpleEvaluator::evaluate_aux(q->left);
//        auto rightGraph = SimpleEvaluator::evaluate_aux(q->right);
//
//        // join left with right
//        return SimpleEvaluator::join(leftGraph, rightGraph);
//
//    }
//
//    return nullptr;
//}

std::vector<RPQTree*> SimpleEvaluator::getLeaves(RPQTree *query) {
    if (query->isLeaf()) {
        return {query};
    }

    std::vector<RPQTree*> result;
    if (query->left) {
        auto rec = getLeaves(query->left);
        result.insert(result.end(), rec.begin(), rec.end());
    }
    if (query->right) {
        auto rec = getLeaves(query->right);
        result.insert(result.end(), rec.begin(), rec.end());
    }

    return result;
}

RPQTree* SimpleEvaluator::optimizeQuery(RPQTree *query) {
    std::vector<RPQTree*> leaves = getLeaves(query);

    while (leaves.size() > 1) {
        uint32_t bestScore = 0;
        RPQTree *bestTree = nullptr;
        int index = -1;

        for (int i = 0; i < leaves.size()-1; ++i) {
            std::string data("/");
            auto *currentTree = new RPQTree(data, leaves[i], leaves[i+1]);
            uint32_t currentScore = est->estimate(currentTree).noPaths;

            if (bestScore == 0 || bestScore > currentScore) {
                bestScore = currentScore;
                bestTree = currentTree;
                index = i;
            }
        }

        leaves.erase(leaves.begin() + index + 1);
        leaves[index] = bestTree;
    }

    return leaves[0];
}

//cardStatstd::shared_ptr<SimpleGraph> SimpleEvaluator::evaluate_aux(RPQTree *q) {
//
//    // evaluate according to the AST bottom-up
//
//    if(q->isLeaf()) {
//        // project out the label in the AST
//        std::regex directLabel (R"((\d+)\+)");
//        std::regex inverseLabel (R"((\d+)\-)");
//
//        std::smatch matches;
//
//        uint32_t label;
//        bool inverse;
//
//        if(std::regex_search(q->data, matches, directLabel)) {
//            label = (uint32_t) std::stoul(matches[1]);
//            inverse = false;
//        } else if(std::regex_search(q->data, matches, inverseLabel)) {
//            label = (uint32_t) std::stoul(matches[1]);
//            inverse = true;
//        } else {
//            std::cerr << "Label parsing failed!" << std::endl;
//            return nullptr;
//        }
//
//        //return SimpleEvaluator::project(label, inverse, graph);
//        //return SimpleEvaluator::project_agg_index(label, inverse, graph);
//        return SimpleEvaluator::project_exh_index(label, inverse, graph);
//    }
//
//    if(q->isConcat()) {
//
//        // evaluate the children
//        auto leftGraph = SimpleEvaluator::evaluate_aux(q->left);
//        auto rightGraph = SimpleEvaluator::evaluate_aux(q->right);
//
//        // join left with right
//        return SimpleEvaluator::join(leftGraph, rightGraph);
//
//    }
//
//    return nullptr;
//}

std::vector<std::string> SimpleEvaluator::treeToString(RPQTree *q) {
    std::vector<std::string> vec;
    SimpleEvaluator::treeToString(q, vec);
    return vec;
}

void SimpleEvaluator::treeToString(RPQTree *q, std::vector<std::string> &vec) {
    if (q->isLeaf()) {
        vec.push_back(q->data);
    } else {
        SimpleEvaluator::treeToString(q->left, vec);
        SimpleEvaluator::treeToString(q->right, vec);
    }
}

cardStat SimpleEvaluator::evaluate(RPQTree *query) {

    vector <string> paths;
    vector <shared_ptr<SimpleGraph>> projections;
    shared_ptr<SimpleGraph> result = nullptr;

    cout << endl;
    // Initalize a vector with the labels
    paths = SimpleEvaluator::treeToString(query);
//    for (int i=0; i < paths.size(); i++) {
//        cout << paths[i] << " | ";
//    }
//    cout << endl;

    // Project all the labels
    for (int i=0; i < paths.size(); i++) {
        uint32_t label = (uint32_t) std::stoul(paths[i].substr(0, paths[i].length()-1));
        bool inverse = paths[i].at(1) == '-';
        projections.push_back(project(label, inverse, graph));
    }

    while (paths.size() > 1) {

        // Find the cheapest join
        vector <int> estimate;
        for (int i=0; i < paths.size()-1; i++) {
            string path = paths[i] + "/" + paths[i+1];
            auto ea = est->estimate(RPQTree::strToTree(path));
            estimate.push_back(ea.noPaths);
        }

        int minPos = 0;
        for (unsigned i = 0; i < estimate.size(); ++i )
        {
            if (estimate[i] < estimate[minPos]) {
                minPos = i;
            }
        }

        auto merged_leafs = join(projections[minPos], projections[minPos+1]);
        paths.insert(paths.begin() + minPos, paths[minPos] + "/" + paths[minPos+1]);
        paths.erase(paths.begin() + minPos + 1);
        paths.erase(paths.begin() + minPos + 1);
        projections.insert(projections.begin() + minPos, merged_leafs);
        projections.erase(projections.begin() + minPos + 1);
        projections.erase(projections.begin() + minPos + 1);

    }
    return computeStats(projections[0]);

}
