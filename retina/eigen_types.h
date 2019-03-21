//
// Created by geoff on 2019-03-20.
//

#ifndef RETINA_TYPES_H
#define RETINA_TYPES_H

#include <Eigen/Dense>


namespace Eigen{
    // boolean dynamic-size matrix
    typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;
}


#endif //RETINA_TYPES_H
