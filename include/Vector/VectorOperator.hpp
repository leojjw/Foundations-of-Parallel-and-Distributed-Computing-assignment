#ifndef _VECTOROPERATOR_
#define _VECTOROPERATOR_

enum Operations { ADD, MINUS, MUL, DIV, ASIG };

template <class Type, Operations op>
void VectorOperation(Type* a, Type* b, Type* result, unsigned int L) {
    if constexpr (op == ADD) {
        for (int i = 0; i < L; i++)
            result[i] = a[i] + b[i];
    }
    if constexpr (op == MINUS) {
        for (int i = 0; i < L; i++)
            result[i] = a[i] - b[i];
    }
    if constexpr (op == MUL) {
        for (int i = 0; i < L; i++)
            result[i] = a[i] * b[i];
    }
    if constexpr (op == DIV) {
        for (int i = 0; i < L; i++)
            result[i] = a[i] / b[i];
    }
    if constexpr (op == ASIG) {
        for (int i = 0; i < L; i++)
            result[i] = a[i];
    }
}

template <class Type>
void VectorPlusKVector(Type* a, Type* b, Type* result, Type k, unsigned int L) {
    for (int i = 0; i < L; i++)
        result[i] = a[i] + k * b[i];
}

template <class Type>
Type VectorNorm(Type* a, unsigned int L) {
    Type sum = 0;

    for (int i = 0; i < L; i++)
        sum += a[i] * a[i];

    return sum;
}

template <class Type>
Type VectorDotVector(Type* a, Type* b, unsigned int L) {
    Type sum = 0;

    for (int i = 0; i < L; i++)
        sum += a[i] * b[i];

    return sum;
}
#endif