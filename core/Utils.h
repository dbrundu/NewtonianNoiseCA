/*
 *  Created on: 25/02/2022
 *      Author: Davide Brundu
 */

#ifndef MYUTILS_H_
#define MYUTILS_H_


template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


#endif
