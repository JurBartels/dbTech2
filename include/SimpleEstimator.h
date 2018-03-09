//
// Created by Nikolay Yakovets on 2018-02-01.
//

#ifndef QS_SIMPLEESTIMATOR_H
#define QS_SIMPLEESTIMATOR_H

#include "Estimator.h"
#include "SimpleGraph.h"

class SimpleEstimator : public Estimator {

    std::shared_ptr<SimpleGraph> graph;

    std::vector<double> label_count;
    std::vector<std::pair<uint32_t ,uint32_t >> start_end_set_counts;

public:
    explicit SimpleEstimator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEstimator() = default;

    void reduceQuery(RPQTree *q, std::vector<std::pair<uint32_t, bool>> &parsedQuery);
    std::vector<std::pair<uint32_t , uint32_t >> compute_in_end_counts();

    void prepare() override ;
    cardStat estimate(RPQTree *q) override ;

};


#endif //QS_SIMPLEESTIMATOR_H
