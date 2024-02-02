#ifndef __ARRARALLOCATE_
#define __ARRARALLOCATE_

template <class T>
T** Allocate2D(int M, int N) {
    T** tmp;
    tmp = new T*[M];
    for (int i = 0; i < M; i++)
        tmp[i] = new T[N];

    return tmp;
}

template <class T>
void Delte2D(T** buf, int M) {
    for (int i = 0; i < M; i++)
        delete[] buf[i];

    delete[] buf;
}

template <class T>
T* Allocate1D(int M) {
    T* tmp = new T[M];
    return tmp;
}

template <class T>
void Delte1D(T* buf) {
    delete[] buf;
}

#endif