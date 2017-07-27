    //Loop splitting + avoid branching
    void myfunc( float* a, float* b, float* c, size_t n) {

        if (n > 100) {

            for (size_t i=0; i < n/2; ++i) {

                c[i] = a[i] * b[i/2];

            }

            for (size_t i=n/2; i < n; ++i) {

                c[i] = c[i-n/2] + a[i] * b[(i - n/2) / 2 + 1];

            }

        }

    }



