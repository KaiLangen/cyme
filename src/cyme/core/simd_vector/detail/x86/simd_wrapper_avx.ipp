/*
 * CYME, License
 * 
 * Timothee Ewart - Swiss Federal Institute of technology in Lausanne 
 * 
 * Permission is hereby granted, free of charge, to any person or organization
 * obtaining a copy of the software and accompanying documentation covered by
 * this license (the "Software") to use, reproduce, display, distribute,
 * execute, and transmit the Software, and to prepare derivative works of the
 * Software, and to permit third-parties to whom the Software is furnished to
 * do so, all subject to the following:
 * 
 * The copyright notices in the Software and this entire statement, including
 * the above license grant, this restriction and the following disclaimer,
 * must be included in all copies of the Software, in whole or in part, and
 * all derivative works of the Software, unless such copies or derivative
 * works are solely in the form of machine-executable object code generated by
 * a source language processor.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef CYME_SIMD_WRAPPER_AVX_HPP
#define CYME_SIMD_WRAPPER_AVX_HPP

namespace numeric{
    template<>
    inline  simd_trait<int,memory::avx>::register_type _mm_load1<int,memory::avx>( simd_trait<int,memory::avx>::register_type xmm0, const  simd_trait<int,memory::avx>::value_type a){
        return xmm0;
    }

    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_load1<double,memory::avx>( simd_trait<double,memory::avx>::register_type xmm0,  simd_trait<double,memory::avx>::value_type a){
        return _mm256_broadcast_sd(&a);
    }
   
    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_load<double,memory::avx>( simd_trait<double,memory::avx>::register_type xmm0,  simd_trait<double,memory::avx>::const_pointer a){
        return _mm256_load_pd(a);
    }

    template<>
    void _mm_store<double,memory::avx>( simd_trait<double,memory::avx>::register_type xmm0,  simd_trait<double,memory::avx>::pointer a){
        _mm256_store_pd(a,xmm0); 
    }
   
    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_mul<double,memory::avx>( simd_trait<double,memory::avx>::register_type xmm0,  simd_trait<double,memory::avx>::register_type xmm1){
        return _mm256_mul_pd(xmm0, xmm1);
    }
   
    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_div<double,memory::avx>( simd_trait<double,memory::avx>::register_type xmm0,  simd_trait<double,memory::avx>::register_type xmm1){
        return _mm256_div_pd(xmm0, xmm1);
    }
   
    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_add<double,memory::avx>( simd_trait<double,memory::avx>::register_type xmm0,  simd_trait<double,memory::avx>::register_type xmm1){
        return _mm256_add_pd(xmm0, xmm1);
    }

    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_sub<double,memory::avx>( simd_trait<double,memory::avx>::register_type xmm0,  simd_trait<double,memory::avx>::register_type xmm1){
        return _mm256_sub_pd(xmm0, xmm1);
    }

    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_rec<double,memory::avx>(simd_trait<double,memory::avx>::register_type xmm0){
        return _mm256_cvtps_pd(_mm_rcp_ps(_mm256_cvtpd_ps(xmm0))); // 256d --(cast)--> 128s --(cast)--> 256d
    }

    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_neg<double,memory::avx>(simd_trait<double,memory::avx>::register_type xmm0){
        simd_trait<double,memory::avx>::register_type mask =  _mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000));
        return _mm256_xor_pd(xmm0, mask);
    }

    template<>
    inline  simd_trait<int,memory::avx>::register_type _mm_floor<double,memory::avx>(simd_trait<double,memory::avx>::register_type xmm0){
        return _mm256_castsi128_si256(_mm256_cvttpd_epi32(_mm256_floor_pd(xmm0)));
    }

    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_cast<double,memory::avx>(simd_trait<int,memory::avx>::register_type xmm0){
        return  _mm256_cvtepi32_pd(_mm256_castsi256_si128(xmm0));
    }

    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_twok<double,memory::avx>(simd_trait<int,memory::avx>::register_type xmm0){
        // PLEASE TUNE ME
        // ((int + 127) << 23) <=> int to float
        __m128i imm0 =  _mm_shuffle_epi32(_mm_slli_epi32(_mm_add_epi32(_mm256_castsi256_si128(xmm0), _mm_set1_epi32(1023)), 20), _MM_SHUFFLE(1,3,0,2));
        __m128i imm1 =  _mm_slli_epi64(imm0,32);
        imm0 =  _mm_srli_epi64(imm0,32); //mask will be slower because mov + broadcast + and, I need to mask 6 instructions
        imm0 =  _mm_slli_epi64(imm0,32);
        xmm0 =   _mm256_insertf128_si256(xmm0, imm0,0);
        xmm0 =   _mm256_insertf128_si256(xmm0, imm1,1);
        return  _mm256_castsi256_pd(xmm0);
    }

    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_min<double,memory::avx>(simd_trait<double,memory::avx>::register_type xmm0, simd_trait<double,memory::avx>::register_type xmm1){
        return _mm256_min_pd(xmm0,xmm1);
    }

    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_max<double,memory::avx>(simd_trait<double,memory::avx>::register_type xmm0, simd_trait<double,memory::avx>::register_type xmm1){
        return _mm256_max_pd(xmm0,xmm1);
    }

#ifdef __INTEL_COMPILER
    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_exp<double,memory::avx>( simd_trait<double,memory::avx>::register_type xmm0){
        return _mm256_exp_pd(xmm0);
    }

    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_log<double,memory::avx>( simd_trait<double,memory::avx>::register_type xmm0){
        return _mm256_log_pd(xmm0);
    }
#endif
    
#ifdef __FMA__
   template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_fma<double,memory::avx>(simd_trait<double,memory::avx>::register_type xmm0,
                                                                                      simd_trait<double,memory::avx>::register_type xmm1,
                                                                                      simd_trait<double,memory::avx>::register_type xmm2){
        return _mm256_fmadd_pd(xmm0, xmm1, xmm2);
    }

    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_nfma<double,memory::avx>(simd_trait<double,memory::avx>::register_type xmm0,
                                                                                       simd_trait<double,memory::avx>::register_type xmm1,
                                                                                       simd_trait<double,memory::avx>::register_type xmm2){
        return _mm256_fnmadd_pd(xmm0, xmm1, xmm2);
    }

    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_fms<double,memory::avx>(simd_trait<double,memory::avx>::register_type xmm0,
                                                                                      simd_trait<double,memory::avx>::register_type xmm1,
                                                                                      simd_trait<double,memory::avx>::register_type xmm2){
        return _mm256_fmsub_pd(xmm0, xmm1, xmm2);
    }

    template<>
    inline  simd_trait<double,memory::avx>::register_type _mm_nfms<double,memory::avx>(simd_trait<double,memory::avx>::register_type xmm0,
                                                                                       simd_trait<double,memory::avx>::register_type xmm1,
                                                                                       simd_trait<double,memory::avx>::register_type xmm2){
        return _mm256_fnmsub_pd(xmm0, xmm1, xmm2);
    }
#endif //end FMA
    
    template<>
     simd_trait<float,memory::avx>::register_type _mm_load1<float,memory::avx>( simd_trait<float,memory::avx>::register_type xmm0,  simd_trait<float,memory::avx>::value_type a){
        return _mm256_broadcast_ss(&a);
    }
   
    template<>
     simd_trait<float,memory::avx>::register_type _mm_load<float,memory::avx>( simd_trait<float,memory::avx>::register_type xmm0,  simd_trait<float,memory::avx>    ::const_pointer a){
        return _mm256_load_ps(a);
    }

    template<>
    void _mm_store<float,memory::avx>( simd_trait<float,memory::avx>::register_type xmm0,  simd_trait<float,memory::avx>::pointer a){
        _mm256_store_ps(a,xmm0); 
    }
   
    template<>
    inline simd_trait<float,memory::avx>::register_type _mm_mul<float,memory::avx>( simd_trait<float,memory::avx>::register_type xmm0,  simd_trait<float,memory::avx>::register_type xmm1){
        return _mm256_mul_ps(xmm0, xmm1);
    }
   
    template<>
    inline simd_trait<float,memory::avx>::register_type _mm_div<float,memory::avx>( simd_trait<float,memory::avx>::register_type xmm0,  simd_trait<float,memory::avx>::register_type xmm1){
        return _mm256_div_ps(xmm0, xmm1);
    }
   
    template<>
    inline simd_trait<float,memory::avx>::register_type _mm_add<float,memory::avx>( simd_trait<float,memory::avx>::register_type xmm0,  simd_trait<float,memory::avx>::register_type xmm1){
        return _mm256_add_ps(xmm0, xmm1);
    }

    template<>
    inline simd_trait<float,memory::avx>::register_type _mm_sub<float,memory::avx>( simd_trait<float,memory::avx>::register_type xmm0,  simd_trait<float,memory::avx>::register_type xmm1){
        return _mm256_sub_ps(xmm0, xmm1);
    }

    template<>
    inline simd_trait<float,memory::avx>::register_type _mm_rec<float,memory::avx>(simd_trait<float,memory::avx>::register_type xmm0){
        return _mm256_rcp_ps(xmm0);
    }

    template<>
    inline  simd_trait<float,memory::avx>::register_type _mm_neg<float,memory::avx>(simd_trait<float,memory::avx>::register_type xmm0){
        simd_trait<float,memory::avx>::register_type mask =  _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));
        return _mm256_xor_ps(xmm0, mask);
    }

    template<>
    inline  simd_trait<int,memory::avx>::register_type _mm_floor<float,memory::avx>(simd_trait<float,memory::avx>::register_type xmm0){
        return _mm256_cvttps_epi32(_mm256_floor_ps(xmm0));
    }

    template<>
    inline  simd_trait<float,memory::avx>::register_type _mm_cast<float,memory::avx>(simd_trait<int,memory::avx>::register_type xmm0){
        return  _mm256_cvtepi32_ps(xmm0);
    }

    template<>
    inline  simd_trait<float,memory::avx>::register_type _mm_twok<float,memory::avx>(simd_trait<int,memory::avx>::register_type xmm0){
        // ((int + 127) << 23) <=> int to float 
        /* AVX2, return _mm256_castsi256_ps(_mm256_slli_epi32(_mm256_add_epi32(xmm0, _mm256_set1_epi32(127)), 23)); */
        xmm0 = _mm256_insertf128_si256(xmm0, _mm_slli_epi32(_mm_add_epi32(_mm256_extractf128_si256(xmm0,0), _mm_set1_epi32(127)), 23),0);
        xmm0 = _mm256_insertf128_si256(xmm0, _mm_slli_epi32(_mm_add_epi32(_mm256_extractf128_si256(xmm0,1), _mm_set1_epi32(127)), 23),1);
        return  _mm256_castsi256_ps(xmm0);
    }

    template<>
    inline  simd_trait<float,memory::avx>::register_type _mm_min<float,memory::avx>(simd_trait<float,memory::avx>::register_type xmm0, simd_trait<float,memory::avx>::register_type xmm1){
        return _mm256_min_ps(xmm0,xmm1);
    }

    template<>
    inline  simd_trait<float,memory::avx>::register_type _mm_max<float,memory::avx>(simd_trait<float,memory::avx>::register_type xmm0, simd_trait<float,memory::avx>::register_type xmm1){
        return _mm256_max_ps(xmm0,xmm1);
    }

#ifdef  __INTEL_COMPILER
    template<>
    inline simd_trait<float,memory::avx>::register_type _mm_exp<float,memory::avx>( simd_trait<float,memory::avx>::register_type xmm0){
        return _mm256_exp_ps(xmm0);
    }

    template<>
    inline simd_trait<float,memory::avx>::register_type _mm_log<float,memory::avx>( simd_trait<float,memory::avx>::register_type xmm0){
        return _mm256_log_ps(xmm0);
    }
#endif

    
#ifdef __FMA__
    template<>
    inline simd_trait<float,memory::avx>::register_type _mm_fma<float,memory::avx>(simd_trait<float,memory::avx>::register_type xmm0,
                                                                                   simd_trait<float,memory::avx>::register_type xmm1,
                                                                                   simd_trait<float,memory::avx>::register_type xmm2){
        return _mm256_fmadd_ps(xmm0, xmm1, xmm2);
    }

    template<>
    inline simd_trait<float,memory::avx>::register_type _mm_nfma<float,memory::avx>(simd_trait<float,memory::avx>::register_type xmm0,
                                                                                    simd_trait<float,memory::avx>::register_type xmm1,
                                                                                    simd_trait<float,memory::avx>::register_type xmm2){
        return _mm256_fnmadd_ps(xmm0, xmm1, xmm2);
    }

    template<>
    inline simd_trait<float,memory::avx>::register_type _mm_fms<float,memory::avx>(simd_trait<float,memory::avx>::register_type xmm0,
                                                                                   simd_trait<float,memory::avx>::register_type xmm1,
                                                                                   simd_trait<float,memory::avx>::register_type xmm2){
        return _mm256_fmsub_ps(xmm0, xmm1, xmm2);
    }

    template<>
    inline simd_trait<float,memory::avx>::register_type _mm_nfms<float,memory::avx>(simd_trait<float,memory::avx>::register_type xmm0,
                                                                                    simd_trait<float,memory::avx>::register_type xmm1,
                                                                                    simd_trait<float,memory::avx>::register_type xmm2){
        return _mm256_fnmsub_ps(xmm0, xmm1, xmm2);
    }
#endif //end FMA
} //end namespace 

#endif
