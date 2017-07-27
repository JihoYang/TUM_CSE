// IMPORTANT: This file requires C++0x or C++11 support!
#include <iostream>
#include <array>
#include <cstring>

struct Memory {
	/*
	 * Question: How come data is located sequentially in memory?
	 * Answer: No guarantee with classes that they will be sequential in memory, with struct it always(?) does.
	 */
	/*
	 * Question: Does the order of the variables in the class matter in terms of how they will be arranged in memory?
	 * Answer: Yes.
	 */

	char buffer[10];
	const std::array<char,72> readonly;

    Memory() 
    :   readonly( // Member initialization list
        {
         // NOTE: There is a special character in front of the actual message.
         // TODO: What is it's purpose with respect to string operations?
         0x00,0x43,0x6F,0x6E,0x67,0x72,0x61,0x74,0x75,0x6C,
         0x61,0x74,0x69,0x6F,0x6E,0x73,0x2C,0x20,0x79,0x6F,
         0x75,0x20,0x68,0x61,0x76,0x65,0x20,0x73,0x75,0x63,
         0x63,0x65,0x73,0x73,0x66,0x75,0x6C,0x6C,0x79,0x20,
         0x63,0x68,0x61,0x6E,0x67,0x65,0x64,0x20,0x61,0x20,
         0x70,0x72,0x6F,0x67,0x72,0x61,0x6D,0x27,0x73,0x20,
         0x62,0x65,0x68,0x61,0x76,0x69,0x6F,0x75,0x72,0x2E,
         0x0A,0x0
       })
    {
        memset(buffer,0,sizeof(buffer));
    }
};

int main(int argc, char **argv) {
    Memory memory;
 
    // --------------------------------------------------------------------------------------
    // If you run the program for the first time you will notice that
    // no output is shown though there it is clearly not empty.
    // TODO: Increase the variable <offset> until the actual message is shown:
    //
    size_t offset = 1;
    /*
     * Solution:
     * size_t offset = 1;
     */

    std::cout << "readonly content: \"" << memory.readonly.data()+offset << "\"" << std::endl;
    
    // --------------------------------------------------------------------------------------
    std::cout << "buffer content before any changes: \"" << memory.buffer << "\"" << std::endl;
    // --------------------------------------------------------------------------------------

    // --------------------------------------------------------------------------------------
    // TODO: Change the number of bytes given by the variable <overwrite> 
    // until you get the message above plus a bunch of 'A's in front of it.
    size_t overwrite = 11;
    /*
     * Solution:
     * size_t overwrite = 11;
     */

    memset(memory.buffer, 'A', overwrite);
    /*
     * Question: Since 'readonly' is const, why does this overwrite not throw an error?
     * Answer: The const marker is on the variable 'readonly'. We are not accessing the memory
     * 		   of 'readonly' through 'readonly' itself. Since the memory is on the stack,
     * 		   there is nothing to protect it from being overwritten.
     */

    std::cout << "buffer content after some changes: \"" << memory.buffer << "\"" << std::endl;

    // --------------------------------------------------------------------------------------
    return 0;
}
