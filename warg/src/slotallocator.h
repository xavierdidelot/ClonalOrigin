/*
 * slotallocator.h
 *
 *  Created on: Jul 21, 2009
 *      Author: koadman
 */

#ifndef SLOTALLOCATOR_H
#define SLOTALLOCATOR_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <vector>
#include <list>
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <cmath>

namespace weakarg {


/** When more space is needed to store a datatype, the memory pool will grow by this factor */
const double POOL_GROWTH_RATE = 1.6;

/**
 * This class allocates memory according to the slot allocation scheme for
 * fixed size objects.  Each time all slots are full it allocates twice the
 * previous allocation.  If it is unable to allocate twice the previous
 * allocation, it does a binary 'search' for the largest amount of memory it
 * can allocate.
 * The current implementation does not allow memory to
 * be freed once allocated.
 */
template< class T >
class SlotAllocator {
public:
	static SlotAllocator<T>& GetSlotAllocator();
	T* Allocate();
	void Free( T* t );
	void Free( std::vector<T*>& chunk );
	~SlotAllocator(){
		Purge();
	};
	void Purge(){
		std::vector<T*>& data = this->data;
		unsigned& tail_free = this->tail_free;
		unsigned& n_elems = this->n_elems;
		std::vector< T* >& free_list = this->free_list;
		for( unsigned dataI = 0; dataI < data.size(); dataI++ )
			free(data[dataI]);
		data.clear();
		free_list.clear();
		tail_free = 0;
		n_elems = 0;
	}

protected:
	std::vector<T*> data;
	unsigned tail_free;
	unsigned n_elems;	/**< number of T in the most recently allocated block */

	std::vector< T* > free_list;

private:
	SlotAllocator() : tail_free(0), n_elems(0) {};
	SlotAllocator& operator=( SlotAllocator& sa ){ n_elems = sa.n_elems; data = sa.data; tail_free = sa.tail_free; return *this;};
	SlotAllocator( SlotAllocator& sa ){ *this = sa; };

};

template< class T >
inline
SlotAllocator< T >& SlotAllocator< T >::GetSlotAllocator(){
	static SlotAllocator< T >* sa = new SlotAllocator< T >();
	return *sa;
}


template< class T >
inline
T* SlotAllocator< T >::Allocate(){
	T* t_ptr = NULL;

{
	std::vector<T*>& data = this->data;
	unsigned& tail_free = this->tail_free;
	unsigned& n_elems = this->n_elems;
	std::vector< T* >& free_list = this->free_list;
	if( free_list.begin() != free_list.end() ){
		t_ptr = free_list.back();
		free_list.pop_back();
	}else if( tail_free > 0 ){
		int T_index = n_elems - tail_free--;
		t_ptr = &(data.back()[ T_index ]);
	}else{

		// Last resort:
		// increase the size of the data array
		unsigned new_size = round((double)n_elems * POOL_GROWTH_RATE);
		if( new_size == 0 )
			new_size++;
		T* new_data = NULL;
		while( true ){
			try{
				new_data = (T*)malloc(sizeof(T)*new_size);
				break;
			}catch(...){
				new_size = new_size / 2;
				if( new_size == 0 )
					break;
			}
		}
		if( new_data == NULL || new_size == 0 ){
			throw std::out_of_range( "SlotAllocator::Allocate(): Unable to allocate more memory" );
		}
		data.push_back( new_data );
		tail_free = new_size - 1;
		t_ptr = & data.back()[0];
		n_elems = new_size;
	}
}
	return t_ptr;
}

template< class T >
inline
void SlotAllocator< T >::Free( T* t ){
	// for debugging double free
/*	for(size_t i = 0; i < free_list.size(); i++ )
		if( free_list[i] == t )
			std::cerr << "ERROR DOUBLE FREE\n";
*/
	t->~T();
{
	std::vector< T* >& free_list = this->free_list;

	free_list.push_back( t );
}
}

template< class T >
inline
void SlotAllocator< T >::Free( std::vector<T*>& chunk ){
	// for debugging double free
/*	for(size_t i = 0; i < free_list.size(); i++ )
		if( free_list[i] == t )
			std::cerr << "ERROR DOUBLE FREE\n";
*/
	for( size_t i = 0; i < chunk.size(); i++ )
		chunk[i]->~T();
{
	std::vector< T* >& free_list = this->free_list;
	free_list.insert(free_list.end(), chunk.begin(), chunk.end());
}
	chunk.clear();
}

}



#endif /* SLOTALLOCATOR_H */
