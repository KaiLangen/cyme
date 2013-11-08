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
#ifndef CYME_VECTOR_HPP
#define CYME_VECTOR_HPP

#include <vector>
#include "memory/detail/simd.h"
#include "memory/allocator.hpp"
#include "memory/detail/storage.hpp"

namespace memory{
    template<class T, std::size_t M,  memory::order O>
    class block_v{
    };

    /**
     \brief block_v of the memory (partial specialization) is instantiated following  AoS layout
     T is the type, M the size of the subblock_v and N the total number of subblock_v. For AoS  specialization
     datas are contiguous block_v after block_v.b I just instantiate a boost::array of block_v
     */
    template<class T, std::size_t M>
    class block_v<T,M,AoS> : public std::vector<storage<T,M,AoS>, memory::Allocator<storage<T,M,AoS> > > {
    public:
        typedef std::size_t                   size_type;
        typedef T                             value_type;
        typedef value_type&                   reference;
        typedef const value_type&             const_reference;
        typedef storage<T,M,AoS>              storage_type;
        typedef std::vector<storage_type, memory::Allocator<storage_type> >  base_type; //default template seems impossible on partial specialization
        typedef typename  base_type::iterator iterator;

        /**
         \brief Default constructor, the block_v is set up to 0
         */
        explicit block_v(const size_type size, const value_type value)
        :base_type(size,storage_type(value)){
        }

        block_v(block_v<T,M,AoS>const& v):base_type(v.size()){
            std::copy(v.begin(),v.end(),(*this).begin()); // memcopy may be faster ....
        }

        /**
         \brief return the value of the block_v i, element j, write only
         */
        inline reference operator()(size_type i, size_type j){
            BOOST_ASSERT_MSG( i < base_type::size(), "out of range: block_v AoS i" );
            BOOST_ASSERT_MSG( j < M, "out of range: block_v AoS j" );
            return base_type::operator[](i)(j);
        }

        /**
         \brief return the value of the subblock_v i, element j, read only
         */
        inline const_reference operator()(size_type i, size_type j) const{
            BOOST_ASSERT_MSG( i < base_type::size(), "out of range: block_v AoS i" );
            BOOST_ASSERT_MSG( j < M, "out of range: block_v AoS j" );
            return base_type::operator[](i)(j);
        }

        /**
         \brief return the size of basic subblock_v
         */
        static inline size_type size_block() {
            return M;
        }

        /**
        \brief adding a new element at the end of the AoS container, using classical push_back of the mother container
        */
        void push_back(value_type value){
            base_type::push_back(storage_type(value));
        }

        /**
         \brief adding a new element at the beginning of the AoS container, using classical push_back of the mother container
         */
        void push_front(value_type value){
            base_type::insert(base_type::begin(),storage_type(value));
        }
    };

    template<class T, std::size_t M>
    class block_v<T,M,AoSoA> : public std::vector<storage<T,__GETSIMD__()/sizeof(T)*M,AoSoA>, memory::Allocator<storage<T,__GETSIMD__()/sizeof(T)*M,AoSoA> > >{
    public:
        typedef std::size_t                                               size_type;
        typedef T                                                         value_type;
        typedef value_type&                                               reference;
        typedef const value_type&                                         const_reference;
        typedef storage<T,__GETSIMD__()/sizeof(T)*M,AoSoA>                storage_type;
        typedef std::vector<storage_type, memory::Allocator<storage<T,__GETSIMD__()/sizeof(T)*M,AoSoA> > >   base_type;                  //default template seems impossible on partial specialization
        typedef typename  base_type::iterator                             iterator;

        explicit block_v(const size_type size, const value_type value)
        :base_type(size/(__GETSIMD__()/sizeof(T))+1, storage_type(value)),size_cyme(size){
        }

        block_v(block_v<T,M,AoSoA >const& v):base_type(v.size()),size_cyme(v.size()){
            std::copy(v.begin(),v.end(),(*this).begin());
        }

        inline reference operator()(size_type i, size_type j){
           // nothing on i as the original size is destroyed in the constructor 
            BOOST_ASSERT_MSG( i < size(), "out of range: block_v AoS i" );
            BOOST_ASSERT_MSG(     j < M, "out of range: block_v AoSoA j" );
            // Please tune me ! (does it exist an alternative to this ? ^_^
            return base_type::operator[]((i*M+j)/(M*__GETSIMD__()/sizeof(T))) //(i)
            (j*(__GETSIMD__()/sizeof(T)) + i%(__GETSIMD__()/sizeof(T)));      //(j)
        };

        inline const_reference operator()(size_type i, size_type j) const{
           // nothing on i as the original size is destroyed in the constructor 
            BOOST_ASSERT_MSG(     j < M, "out of range: block_v AoSoA j" );
            // Please tune me ! (does it exist an alternative to this ? ^_^
            return base_type::operator[]((i*M+j)/(M*__GETSIMD__()/sizeof(T))) //(i)
            (j*(__GETSIMD__()/sizeof(T)) + i%(__GETSIMD__()/sizeof(T)));      //(j)
        };

        static inline size_type size_block() {
            return M;
        }

        void resize(size_type n){
            return base_type::resize(n/(__GETSIMD__()/sizeof(T))+1);
        }

        void reserve(size_type n){
            return base_type::reserve(n/(__GETSIMD__()/sizeof(T))+1);
        }

        const size_type size() const{
            return size_cyme;
        }


        /**
         \brief adding a new element at the end of the AoSoA container, first check the size if pb increase
         */
        void push_back(value_type value){
            if((this->size_cyme/(__GETSIMD__()/sizeof(T))+1) > base_type::size())
                (*this).resize(size_cyme);

            for(size_type j=0; j<M; ++j)
                (*this)(this->size_cyme,j)=value; // I prefer (*this) than operator()

            ++size_cyme;
        }


        /**
         \brief adding a new element at the beginning of the AoSoA container, first check the size if pb increase, 
          and copy element one by one, very slow ....
         */
        void push_front(value_type value){
            BOOST_ASSERT_MSG(true, " push_front is VERY SLOW for AoSoA container " );
            if((this->size_cyme/(__GETSIMD__()/sizeof(T))+1) > base_type::size())
                (*this).resize(size_cyme);

            for(size_type i=size_cyme; i>0; --i){ // reorder coeff one by one it is slow
                BOOST_ASSERT_MSG( i == 0, " mistake push_front ! Debug ! " );
                for(size_type j=0; j<M; ++j)
                    (*this)(i,j)=(*this)(i-1,j); // the -1 makes the translation
            }

            for(size_type j=0; j<M; ++j)
                (*this)(0,j)=value; // I prefer (*this) than operator()

            ++size_cyme;
        }

    private:
        size_type size_cyme; // it is the same than the size() for AoS
    };
} //end namespace memory

namespace cyme {
    /*
     \brief  This class facilitates the creation of an array of synapses (or whatever), the condition the class
     must encapsulate the basic type (value_type) and the size (value_size) of the basic object under the form:
     template <class T>
     class example{
     typedef T value_type;
     static const int value_size = 5;
     }
     */
    template<class T, memory::order O>
    class vector : public memory::block_v<typename T::value_type,  T::value_size, O>{
    public:
        typedef typename T::value_type value_type;

        explicit vector(const size_t size = 1, const value_type value = value_type())
        :memory::block_v<value_type, T::value_size, O>(size, value){
        }

        vector(vector const& a):memory::block_v<typename T::value_type,  T::value_size, O>(a){
        }
   };
}

#endif
