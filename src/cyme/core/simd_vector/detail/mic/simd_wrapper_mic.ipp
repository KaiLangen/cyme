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

#ifndef CYME_SIMD_WRAPPER_MIC_HPP
#define CYME_SIMD_WRAPPER_MIC_HPP

namespace numeric{

    template<>
    inline  simd_trait<int,memory::mic>::register_type _mm_load1<int,memory::mic>( simd_trait<int,memory::mic>::register_type xmm0, const  simd_trait<int,memory::mic>::value_type a){
        return xmm0;
    }

    /*  -------------------------------------------------------------------------------------------------------------------- double */
    template<>
    inline typename simd_trait<double,memory::mic>::register_type _mm_load1<double,memory::mic>(typename simd_trait<double,memory::mic>::register_type xmm0, typename simd_trait<double,memory::mic>::value_type a){
        return _mm512_set1_pd(a);
    }
   
    template<>
    inline typename simd_trait<double,memory::mic>::register_type _mm_load<double,memory::mic>(typename simd_trait<double,memory::mic>::register_type xmm0, typename simd_trait<double,memory::mic>::const_pointer a){
        return _mm512_load_pd(a);
    }

    template<>
    void _mm_store<double,memory::mic>(typename simd_trait<double,memory::mic>::register_type xmm0, typename simd_trait<double,memory::mic>::pointer a){
        _mm512_store_pd(a,xmm0); 
    }
   
    template<>
    inline typename simd_trait<double,memory::mic>::register_type _mm_mul<double,memory::mic>(typename simd_trait<double,memory::mic>::register_type xmm0, typename simd_trait<double,memory::mic>::register_type xmm1){
        return _mm512_mul_pd(xmm0, xmm1);
    }
   
    template<>
    inline typename simd_trait<double,memory::mic>::register_type _mm_div<double,memory::mic>(typename simd_trait<double,memory::mic>::register_type xmm0, typename simd_trait<double,memory::mic>::register_type xmm1){
        return _mm512_div_pd(xmm0, xmm1);
    }
   
    template<>
    inline typename simd_trait<double,memory::mic>::register_type _mm_add<double,memory::mic>(typename simd_trait<double,memory::mic>::register_type xmm0, typename simd_trait<double,memory::mic>::register_type xmm1){
        return _mm512_add_pd(xmm0, xmm1);
    }

    template<>
    inline typename simd_trait<double,memory::mic>::register_type _mm_sub<double,memory::mic>(typename simd_trait<double,memory::mic>::register_type xmm0, typename simd_trait<double,memory::mic>::register_type xmm1){
        return _mm512_sub_pd(xmm0, xmm1);
    }

    template<>
    inline  simd_trait<double,memory::mic>::register_type _mm_neg<double,memory::mic>(simd_trait<double,memory::mic>::register_type xmm0){
        return (xmm0);
    };

    template<>
    inline  simd_trait<int,memory::mic>::register_type _mm_floor<double,memory::mic>(simd_trait<double,memory::mic>::register_type xmm0){
        return (xmm0);
    };

    template<>
    inline  simd_trait<double,memory::mic>::register_type _mm_cast<double,memory::mic>(simd_trait<int,memory::mic>::register_type xmm0){
        return  (xmm0);
    };

    template<>
    inline  simd_trait<double,memory::mic>::register_type _mm_twok<double,memory::mic>(simd_trait<int,memory::mic>::register_type xmm0){
        return  (xmm0);
    };

    template<>
    inline  simd_trait<double,memory::mic>::register_type _mm_min<double,memory::mic>(simd_trait<double,memory::mic>::register_type xmm0, simd_trait<double,memory::mic>::register_type xmm1){
        return _mm512_min_pd(xmm0,xmm1);
    };

    template<>
    inline  simd_trait<double,memory::mic>::register_type _mm_max<double,memory::mic>(simd_trait<double,memory::mic>::register_type xmm0, simd_trait<double,memory::mic>::register_type xmm1){
        return _mm512_max_pd(xmm0,xmm1);
    };

    template<>
    inline typename simd_trait<double,memory::mic>::register_type _mm_exp<double,memory::mic>(typename simd_trait<double,memory::mic>::register_type xmm0){
        return _mm512_exp_pd(xmm0);
    }

    template<>
    inline typename simd_trait<double,memory::mic>::register_type _mm_log<double,memory::mic>(typename simd_trait<double,memory::mic>::register_type xmm0){
        return _mm512_log_pd(xmm0);
    }

    template<>
    inline simd_trait<double,memory::mic>::register_type _mm_rec<double,memory::mic>(simd_trait<double,memory::mic>::register_type xmm0){
        return (xmm0);
    };

   template<>
    inline typename simd_trait<double,memory::mic>::register_type _mm_fma<double,memory::mic>(typename simd_trait<double,memory::mic>::register_type xmm0, typename simd_trait<double,memory::mic>::register_type xmm1, typename simd_trait<double,memory::mic>::register_type xmm2){
        return _mm512_fmadd_pd(xmm0, xmm1, xmm2);
    }

    template<>
    inline typename simd_trait<double,memory::mic>::register_type _mm_fms<double,memory::mic>(typename simd_trait<double,memory::mic>::register_type xmm0, typename simd_trait<double,memory::mic>::register_type xmm1, typename simd_trait<double,memory::mic>::register_type xmm2){
        return _mm512_fmsub_pd(xmm0, xmm1, xmm2);
    }

    template<>
    inline typename simd_trait<double,memory::mic>::register_type _mm_nfms<double,memory::mic>(typename simd_trait<double,memory::mic>::register_type xmm0, typename simd_trait<double,memory::mic>::register_type xmm1, typename simd_trait<double,memory::mic>::register_type xmm2){
        return _mm512_fnmsub_pd(xmm0, xmm1, xmm2);
    }

    template<>
    inline typename simd_trait<double,memory::mic>::register_type _mm_nfma<double,memory::mic>(typename simd_trait<double,memory::mic>::register_type xmm0, typename simd_trait<double,memory::mic>::register_type xmm1, typename simd_trait<double,memory::mic>::register_type xmm2){
        return _mm512_fnmadd_pd(xmm0, xmm1, xmm2);
    }

    /*  -------------------------------------------------------------------------------------------------------------------- float */
    template<>
    typename simd_trait<float,memory::mic>::register_type _mm_load1<float,memory::mic>(typename simd_trait<float,memory::mic>::register_type xmm0, typename simd_trait<float,memory::mic>::value_type a){
        return _mm512_set1_ps(a);
    }
   
    template<>
    typename simd_trait<float,memory::mic>::register_type _mm_load<float,memory::mic>(typename simd_trait<float,memory::mic>::register_type xmm0, typename simd_trait<float,memory::mic>    ::const_pointer a){
        return _mm512_load_ps(a);
    }

    template<>
    void _mm_store<float,memory::mic>(typename simd_trait<float,memory::mic>::register_type xmm0, typename simd_trait<float,memory::mic>::pointer a){
        _mm512_store_ps(a,xmm0); 
    }
   
    template<>
    inline typename simd_trait<float,memory::mic>::register_type _mm_mul<float,memory::mic>(typename simd_trait<float,memory::mic>::register_type xmm0, typename simd_trait<float,memory::mic>::register_type xmm1){
        return _mm512_mul_ps(xmm0, xmm1);
    }
   
    template<>
    inline typename simd_trait<float,memory::mic>::register_type _mm_div<float,memory::mic>(typename simd_trait<float,memory::mic>::register_type xmm0, typename simd_trait<float,memory::mic>::register_type xmm1){
        return _mm512_div_ps(xmm0, xmm1);
    }
   
    template<>
    inline typename simd_trait<float,memory::mic>::register_type _mm_add<float,memory::mic>(typename simd_trait<float,memory::mic>::register_type xmm0, typename simd_trait<float,memory::mic>::register_type xmm1){
        return _mm512_add_ps(xmm0, xmm1);
    }

    template<>
    inline typename simd_trait<float,memory::mic>::register_type _mm_sub<float,memory::mic>(typename simd_trait<float,memory::mic>::register_type xmm0, typename simd_trait<float,memory::mic>::register_type xmm1){
        return _mm512_sub_ps(xmm0, xmm1);
    }

    template<>
    inline  simd_trait<float,memory::mic>::register_type _mm_neg<float,memory::mic>(simd_trait<float,memory::mic>::register_type xmm0){

        return (xmm0);
    };

    template<>
    inline  simd_trait<int,memory::mic>::register_type _mm_floor<float,memory::mic>(simd_trait<float,memory::mic>::register_type xmm0){
        return (xmm0);
    };

    template<>
    inline  simd_trait<float,memory::mic>::register_type _mm_cast<float,memory::mic>(simd_trait<int,memory::mic>::register_type xmm0){
        return  (xmm0);
    };

    template<>
    inline  simd_trait<float,memory::mic>::register_type _mm_twok<float,memory::mic>(simd_trait<int,memory::mic>::register_type xmm0){
        return  (xmm0);
    };

    template<>
    inline  simd_trait<float,memory::mic>::register_type _mm_min<float,memory::mic>(simd_trait<float,memory::mic>::register_type xmm0, simd_trait<float,memory::mic>::register_type xmm1){
        return _mm512_min_ps(xmm0,xmm1);
    };

    template<>
    inline  simd_trait<float,memory::mic>::register_type _mm_max<float,memory::mic>(simd_trait<float,memory::mic>::register_type xmm0, simd_trait<float,memory::mic>::register_type xmm1){
        return _mm512_max_ps(xmm0,xmm1);
    };



    template<>
    inline simd_trait<float,memory::mic>::register_type _mm_rec<float,memory::mic>(simd_trait<float,memory::mic>::register_type xmm0){
        return _mm512_rcp23_ps(xmm0);
    };

    template<>
    inline typename simd_trait<float,memory::mic>::register_type _mm_exp<float,memory::mic>(typename simd_trait<float,memory::mic>::register_type xmm0){
        return _mm512_exp_ps(xmm0);
    }
    
    template<>
    inline typename simd_trait<float,memory::mic>::register_type _mm_fma<float,memory::mic>(typename simd_trait<float,memory::mic>::register_type xmm0, typename simd_trait<float,memory::mic>::register_type xmm1, typename simd_trait<float,memory::mic>::register_type xmm2){
        return _mm512_fmadd_ps(xmm0, xmm1, xmm2);
    }

    template<>
    inline typename simd_trait<float,memory::mic>::register_type _mm_fms<float,memory::mic>(typename simd_trait<float,memory::mic>::register_type xmm0, typename simd_trait<float,memory::mic>::register_type xmm1, typename simd_trait<float,memory::mic>::register_type xmm2){
        return _mm512_fmsub_ps(xmm0, xmm1, xmm2);
    }

    template<>
    inline typename simd_trait<float,memory::mic>::register_type _mm_nfms<float,memory::mic>(typename simd_trait<float,memory::mic>::register_type xmm0, typename simd_trait<float,memory::mic>::register_type xmm1, typename simd_trait<float,memory::mic>::register_type xmm2){
        return _mm512_fnmsub_ps(xmm0, xmm1, xmm2);
    }

    template<>
    inline typename simd_trait<float,memory::mic>::register_type _mm_nfma<float,memory::mic>(typename simd_trait<float,memory::mic>::register_type xmm0, typename simd_trait<float,memory::mic>::register_type xmm1, typename simd_trait<float,memory::mic>::register_type xmm2){
        return _mm512_fnmadd_ps(xmm0, xmm1, xmm2);
    }

} //end namespace 

#endif 
