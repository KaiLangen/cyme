/*
 * Cyme - simd_exp.hpp, Copyright (c), 2014,
 * Timothee Ewart - Swiss Federal Institute of technology in Lausanne,
 * timothee.ewart@epfl.ch,
 * Kai Langen,
 * kai.langen@usask.ca,
 * All rights reserved.
 * This file is part of Cyme <https://github.com/BlueBrain/cyme>
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

/**
 * @file cyme/core/simd_vector/math/simd_gather.ipp
 * Implements gather for vec_simd class
 */

#ifndef CYME_SIMD_GATHER_IPP
#define CYME_SIMD_GATHER_IPP

namespace cyme{

    template<class T, cyme::simd O, int N>
    struct cyme_gather{
        static forceinline vec_simd<T,O,N> gather(const_pointer a, vec_simd<std::size_t,O,N> const& v){
            const size_t size = elems_helper<T,N>::size;
            T elems[size] __attribute__((aligned(static_cast<size_t>(cyme::trait_register<T,cyme::__GETSIMD__()>::size))));
            for(size_t i = 0; i < size; i++){
//		assert(idx[i] < length);
		elems[i] = a[v[i]];
            }	
            xmm = _mm_load<typename simd_trait<T,O,N>::value_type,O,N>(elems);
            return *this;

    } 
    /** Free function for call the vendor gather */
    template<class T,cyme::simd O, int N>
    forceinline vec_simd<T,O,N> gather_v(const_pointer a, vec_simd<std::size_t,O,N> const& v){
        vec_simd<T,O,N> nrv(_mm_gather<T,O,N>(a,v.xmm));
        return nrv;
    }

    /** Function object for the vendor gather algorithm */
    template<class T, cyme::simd O, int N, std::size_t n>
    struct Vendor_gather{
        static forceinline vec_simd<T,O,N> gather(const_pointer a, vec_simd<std::size_t,O,N> const& v){
            return gather_v(a,v); /* call vendor wrapper */
        }
    };

    /** Selector for the gather algorithm (vendor or cyme implementation) */
    template<class T, cyme::simd O, int N, class Solver = cyme_gather<T,O,N,n> >
    struct Selector_gather{
         static forceinline vec_simd<T,O,N> gather(const_pointer a, vec_simd<std::size_t,O,N> const& v){
               x = Solver::gather(a,v);
               return x;
         }
    };

    /** free function for gather */
    template<class T,cyme::simd O, int N>
    forceinline vec_simd<T,O,N> gather(const_pointer a, vec_simd<std::size_t,O,N> const& v){
        return Selector_gather<T,O,N>::gather(a,v);
    }
}
#endif
