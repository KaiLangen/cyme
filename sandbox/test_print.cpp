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

template<class B>
void init(B& b){
    for(std::size_t i=0; i<b.size(); ++i)
        for(std::size_t j=0; j<b.size_block(); ++j){
             b(i,j) = rand()%10;
        };
}

int main(int argc, char *argv[]){
    cyme::vector<synapse<size_t,8>,cyme::AoSoA> b(32);
    cyme::vector<synapse<double,8>,cyme::AoSoA> c(32);

    init(b);
    init(c);
    cyme::vector<synapse<size_t,8>,cyme::AoSoA>::const_iterator it1 = b.begin();
    cyme::vector<synapse<double,8>,cyme::AoSoA>::const_iterator it2 = c.begin();


    std::cout << (*it2)[(*it1)[0]] << std::endl;
//    std::cout << (*it1)[0] << std::endl;
//    std::cout << (*it2)[5] << std::endl;
    return 0;
};
