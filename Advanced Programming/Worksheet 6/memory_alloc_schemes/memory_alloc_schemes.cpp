#include <iostream>
#include <stdlib.h>

#define array_size 100

int global_array[array_size];
void memory_accesses_read(int* pointer, int start, int stop);
void memory_accesses_write(int* pointer, int start, int stop);


int main(int argc, char** argv){

	static int static_array[array_size];
	int stack_array[array_size];
	int* heap_array = new int [array_size];

	//Print out the address locations of the different arrays
	std::cout << "static address: " << reinterpret_cast<void*>(static_array) << std::endl;
	std::cout << "global address: " << reinterpret_cast<void*>(global_array) << std::endl;
	std::cout << "stack address: " << reinterpret_cast<void*>(stack_array) << std::endl;
	std::cout << "heap address: " << reinterpret_cast<void*>(heap_array) << std::endl;

	delete[] heap_array;

	int start = 0;
	int end = array_size;
	int access = 0;
	while(true){
		std::cout << "[start, end, access]=" << start << "," << end << "," << access << std::endl;
		
		//STATIC
		//memory_accesses_write(static_array,start,end,access); 
		//memory_accesses_read(static_array,start,end);
		
		//GLOBAL
		//memory_accesses_write(global_array,start,end,access); 
		//memory_accesses_read(global_array,start,end);
		
		//STACK
		//memory_accesses_write(stack_array,start,end,access);
		//memory_accesses_read(stack_array,start,end);
		
		//Deleted HEAP
		memory_accesses_write(heap_array,start,end);
		memory_accesses_read(heap_array,start,end);

		start--;
		//end++;
		access++;
	}

}

/**
 * Iterates from [start <= i < stop] and reads value in *(pointer+i)
 * @param pointer
 * @param start
 * @param stop
 */
void memory_accesses_read(int* pointer, int start, int stop){
	int counter = 0;
	std::cout << "#Reading...";
	int read;
	for(int i = start; i<=stop; i++){
		read = *(pointer+i);
		//std::cout << counter << "|" << access << std::endl;
		counter++;
	}
	std::cout << counter << " times. Last read: " << read << std::endl;
}

/**
 * Iterates from [start <= i < stop] and writes random value in *(pointer+i)
 * @param pointer
 * @param start
 * @param stop
 * @param access
 */
void memory_accesses_write(int* pointer, int start, int stop){
	std::cout << "#Writing...";
	int counter = 0;
	for(int i = start; i<=stop; i++){
		*(pointer+i)=rand();
		counter++;
	}
	std::cout << counter << " times." << std::endl;
}
