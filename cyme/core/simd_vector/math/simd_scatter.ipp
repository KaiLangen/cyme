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

#ifndef CYME_SIMD_SCATTER_IPP
#define CYME_SIMD_SCATTER_IPP

#define integer_type typename trait_integer<T>::integer

namespace cyme{
    template<class T,cyme::simd O, int N>
    struct cyme_scatter{
	static void scatter(typename simd_trait<T,O,N>::pointer dst, vec_simd<T,O,N> const& src, vec_simd<integer_type,O,N> const& v){
            const std::size_t size = elems_helper<T,N>::size;
            T elems[size] __attribute__((aligned(static_cast<std::size_t>(cyme::trait_register<T,cyme::__GETSIMD__()>::size))));
            integer_type V[size] __attribute__((aligned(static_cast<std::size_t>(cyme::trait_register<T,cyme::__GETSIMD__()>::size))));
	    _mm_store<T,O,N>(src.xmm,elems);
	    _mm_store<integer_type,O,N>(v.xmm,V);
            for(std::size_t i = 0; i < size; i++){
//		assert(V[i] < size);
		dst[V[i]] = elems[i];
            }	
	}
    };
 
    /** Free function for call the vendor scatter */
    template<class T,cyme::simd O, int N>
    void scatter_v(typename simd_trait<T,O,N>::pointer dst, vec_simd<T,O,N> const& src, vec_simd<integer_type,O,N> const& v){
        _mm_scatter<T,O,N>(dst,src.xmm,v.xmm);
    }

    /** Function object for the vendor scatter algorithm */
    template<class T,cyme::simd O, int N>
    struct Vendor_scatter{
	static void scatter(typename simd_trait<T,O,N>::pointer dst, vec_simd<T,O,N> const& src, vec_simd<integer_type,O,N> const& v){
            scatter_v(dst,src,v); /* call vendor wrapper */
        }
    };

    /** Selector for the scatter algorithm (vendor or cyme implementation) */
    template<class T,cyme::simd O, int N, class Solver = cyme_scatter<T,O,N> >
    struct Selector_scatter{
	static void scatter(typename simd_trait<T,O,N>::pointer dst, vec_simd<T,O,N> const& src, vec_simd<integer_type,O,N> const& v){
               Solver::scatter(dst,src,v);
         }
    };

    /** Free function for scatter */
    template<class T,cyme::simd O, int N>
    void scatter(typename simd_trait<T,O,N>::pointer dst, vec_simd<T,O,N> const& src,
							  vec_simd<typename trait_integer<T>::integer,O,N> const& v){
        Selector_scatter<T,O,N>::scatter(dst,src,v);
    }
}

#undef integer_type
#endif
