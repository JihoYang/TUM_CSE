#include <iostream>

const int N = 10;

class AbstractVector {
    public:
        AbstractVector() {
            std::cout << "AbstractVector: constructor called" << std::endl;
        }

        ~AbstractVector() { // ---- TASK 2 ----
            std::cout << "AbstractVector: destructor called" << std::endl;
        }

        void print() {
            std::cout << "--- PRINT: Abstract vector content is abstract" << std::endl;
        }
};

class FixedVector : public AbstractVector {
    public:
        FixedVector() {
            std::cout << "FixedVector: constructor called" << std::endl;
            for (int i=0; i< N; ++i) {
                data[i] = 0;
            }
        }
 
        ~FixedVector() { // ---- TASK 2 ----
            std::cout << "FixedVector: destructor called" << std::endl;
        }

        void copy(FixedVector target) {
            std::cout << "--- copy data around!" << std::endl;
            for (int i=0; i < N ; ++i) {
               target.data[i] = data[i];
            }
        }

        void print() {
            std::cout << "--- PRINT: Fixed Vector content:";
            for (int i=0; i < N; ++i) {
                std::cout << " " << data[i];
            }
            std::cout << std::endl;
        }

    private:
        int data[N];
};

class SequenceVector : public FixedVector {
    public:
        SequenceVector() {
            std::cout << "SequenceVector: constructor called" << std::endl;
 
            // TODO: allow SequenceVector to directly access the data array in FixedVector
            // then fill the array with values from [0..N-1]
        }
 
        ~SequenceVector() { // ---- TASK 2 ----
            std::cout << "SequenceVector: destructor called" << std::endl;
        }
};

void heap_initialization() {
    std::cout << "heap initialization: " << std::endl;
    std::cout << "----" << std::endl;
    {
        AbstractVector *fixedVector = new FixedVector();
        delete fixedVector;
    }
    std::cout << "----" << std::endl;
    {
        AbstractVector *sequenceVector = new SequenceVector();
        delete sequenceVector;
    }
    std::cout << "----" << std::endl;
}

void normal_initialization() {
    std::cout << "normal initialization: " << std::endl;
    std::cout << "----" << std::endl;
    {
        FixedVector fixedVector;
    }
    std::cout << "----" << std::endl;
    {
        SequenceVector seqVector;
    }
    std::cout << "----" << std::endl;
    std::cout << std::endl;
}

int main(int argc, char **argv) {
    std::cout << "========== TASK 1 =============" << std::endl;
    // TODO 1: direct access for derived classes only
    SequenceVector seqVector;
    std::cout << "this should be a sequence: ";
    seqVector.print();
 
    /* ----------------------------------------------------------- */

    std::cout << "========== TASK 2 =============" << std::endl;
    // TODO 2: destructors and inheritance: 
    // heap version looks different? Why? Fix it!
    normal_initialization();
    heap_initialization();
  
    /* ----------------------------------------------------------- */
 
    std::cout << "========== TASK 3 =============" << std::endl;
    // TODO 3: virtual functions matter
    AbstractVector *absVector = &seqVector;
    absVector->print();

    /* ----------------------------------------------------------- */

    std::cout << "========== TASK 4 =============" << std::endl;
    // TODO 4: call by something ?
    FixedVector newVector;
    seqVector.copy(newVector);
    
    std::cout << "expected values: ";
    seqVector.print();
    std::cout << "actual values: ";
    newVector.print();
 
    return 0;
}
