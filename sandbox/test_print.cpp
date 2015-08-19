/*
 * Cyme - test_print.cpp, Copyright (c), 2014,
 * Timothee Ewart - Swiss Federal Institute of technology in Lausanne,
 * timothee.ewart@epfl.ch,
 * All rights reserved.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library.
 */

#include <cstddef>
#include <cyme/cyme.h>

template<class T, size_t M>
struct synapse{
   typedef T value_type;
   static const size_t value_size = M;
};

template<class T>
T random();

template<>
int random<int>(){
    return rand()%6;
};

template<>
float random<float>(){
    return drand48();
};

template<class Ba, class Bb> // m and n are differents into the block that why I passe like argument
void init(Ba& block_a, Bb& block_b){
    typedef typename Ba::value_type value_type;
    for(std::size_t i=0; i<block_a.size(); ++i)
        for(std::size_t j=0; j<block_a.size_block(); ++j){
            block_a(i,j) = random<value_type>();
            block_b(i,j) = block_a(i,j);
        }
}

int main(int argc, char *argv[]){
    cyme::vector<synapse<int,6>,cyme::AoSoA> b_AOSOA(32);
    cyme::vector<synapse<float,6>,cyme::AoSoA> c_AOSOA(32);

    cyme::vector<synapse<int,6>,cyme::AoS> b_AOS(32);
    cyme::vector<synapse<float,6>,cyme::AoS> c_AOS(32);
    // The init was wrong, and I shoud write a special functor for AOSOA
    // we need a AOS to do the initialisation
    // it is super messy and hgard to understand
    init(b_AOS,b_AOSOA);
    init(c_AOS,c_AOSOA);

    {
        cyme::vector<synapse<int,6>,cyme::AoSoA>::const_iterator it1 = b_AOSOA.begin();
        cyme::vector<synapse<float,6>,cyme::AoSoA>::const_iterator it2 = c_AOSOA.begin();

        std::cout << (*it2)[0] << std::endl;
        std::cout << (*it1)[0] << std::endl;
        std::cout << (*it2)[(*it1)[0]] << std::endl;
    }
    {
        cyme::vector<synapse<int,6>,cyme::AoS>::const_iterator it1 = b_AOS.begin();
        cyme::vector<synapse<float,6>,cyme::AoS>::const_iterator it2 = c_AOS.begin();

        std::cout << (*it2)[0] << std::endl;
        std::cout << (*it1)[0] << std::endl;
        std::cout << (*it2)[(*it1)[0]] << std::endl;
    }

    return 0;
};

