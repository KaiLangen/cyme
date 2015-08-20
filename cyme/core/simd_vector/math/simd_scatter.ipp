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
 * @file cyme/core/simd_vector/math/simd_scatter.ipp
 * Implements scatter for vec_simd class
 */

#ifndef CYME_SIMD_GATHER_IPP
#define CYME_SIMD_GATHER_IPP

#define __UINT__ typename trait_integer<T>::integer

namespace cyme{
    template<class T,cyme::simd O, int N>
    struct cyme_scatter{
        static forceinline vec_simd<T,O,N> scatter(typename simd_trait<T,O,N>::const_pointer a, vec_simd<__UINT__,O,N> const& v){
            const std::size_t size = elems_helper<T,N>::size;
            T elems[size] __attribute__((aligned(static_cast<std::size_t>(cyme::trait_register<T,cyme::__GETSIMD__()>::size))));
            __UINT__ V[size] __attribute__((aligned(static_cast<std::size_t>(cyme::trait_register<T,cyme::__GETSIMD__()>::size))));
	    _mm_store<__UINT__,O,N>(v.xmm,V);
            for(std::size_t i = 0; i < size; i++){
		assert(V[i] < size);
		elems[V[i]] = a[i];
            }	
            vec_simd<T,O,N> nrv(_mm_load<typename simd_trait<T,O,N>::value_type,O,N>(elems));
            return nrv;
	}
    };
 
    /** Free function for call the vendor scatter */
    template<class T,cyme::simd O, int N>
    forceinline vec_simd<T,O,N> scatter_v(typename simd_trait<T,O,N>::const_pointer a, vec_simd<__UINT__,O,N> const& v){
        vec_simd<T,O,N> nrv(_mm_scatter<T,O,N>(a,v.xmm));
        return nrv;
    }

    /** Function object for the vendor scatter algorithm */
    template<class T,cyme::simd O, int N>
    struct Vendor_scatter{
        static forceinline vec_simd<T,O,N> scatter(typename simd_trait<T,O,N>::const_pointer a, vec_simd<__UINT__,O,N> const& v){
            return scatter_v(a,v); /* call vendor wrapper */
        }
    };

    /** Selector for the scatter algorithm (vendor or cyme implementation) */
    template<class T,cyme::simd O, int N, class Solver = cyme_scatter<T,O,N> >
    struct Selector_scatter{
         static forceinline vec_simd<T,O,N> scatter(typename simd_trait<T,O,N>::const_pointer a, vec_simd<__UINT__,O,N> const& v){
               return Solver::scatter(a,v);
         }
    };

    /** free function for scatter */
    template<class T,cyme::simd O, int N>
    forceinline vec_simd<T,O,N> scatter(typename simd_trait<T,O,N>::const_pointer a, vec_simd<__UINT__,O,N> const& v){
        return Selector_scatter<T,O,N>::scatter(a,v);
    }
}

#undef __UINT__
#endif
