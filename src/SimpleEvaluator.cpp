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

std::shared_ptr<SimpleGraph> SimpleEvaluator::evaluate_aux(RPQTree *q) {

    // evaluate according to the AST bottom-up

    if(q->isLeaf()) {
        // project out the label in the AST
        std::regex directLabel (R"((\d+)\+)");
        std::regex inverseLabel (R"((\d+)\-)");

        std::smatch matches;

        uint32_t label;
        bool inverse;

        if(std::regex_search(q->data, matches, directLabel)) {
            label = (uint32_t) std::stoul(matches[1]);
            inverse = false;
        } else if(std::regex_search(q->data, matches, inverseLabel)) {
            label = (uint32_t) std::stoul(matches[1]);
            inverse = true;
        } else {
            std::cerr << "Label parsing failed!" << std::endl;
            return nullptr;
        }

        return SimpleEvaluator::project(label, inverse, graph);
    }

    if(q->isConcat()) {

        // evaluate the children
        auto leftGraph = SimpleEvaluator::evaluate_aux(q->left);
        auto rightGraph = SimpleEvaluator::evaluate_aux(q->right);

        // join left with right
        return SimpleEvaluator::join(leftGraph, rightGraph);

    }

    return nullptr;
}

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

cardStat SimpleEvaluator::evaluate(RPQTree *query) {
    vector <string> paths;
    vector <shared_ptr<SimpleGraph>> projections;
    shared_ptr<SimpleGraph> result = nullptr;

    // Initalize a vector with the labels
    paths = SimpleEvaluator::treeToString(query);

    // Project all the labels
    for (int i=0; i < paths.size(); i++) {
        uint32_t label = (uint32_t) std::stoul(paths[i].substr(0, 1));
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

