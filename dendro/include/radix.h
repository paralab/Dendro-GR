//
// Created by milinda on 5/26/16.
//

#ifndef SFCSORTBENCH_RADIX_H
#define SFCSORTBENCH_RADIX_H

// A utility function to get maximum value in arr[]
template <typename T>
inline T getMax(T* arr, int n)
{
    T mx = arr[0];
    for (unsigned int i = 1; i < n; i++)
        if (arr[i] > mx)
            mx = arr[i];
    return mx;
}

// A function to do counting sort of arr[] according to
// the digit represented by exp.
template <typename T>
inline void countSort(T* arr, int n, int exp)
{
    unsigned int output[n]; // output array
    unsigned int i, count[10] = {0};

    // Store count of occurrences in count[]
    for (i = 0; i < n; i++)
        count[ (arr[i]/exp)%10 ]++;

    // Change count[i] so that count[i] now contains actual
    //  position of this digit in output[]
    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];

    // Build the output array
    for (i = n - 1; i >= 0; i--)
    {
        output[count[ (arr[i]/exp)%10 ] - 1] = arr[i];
        count[ (arr[i]/exp)%10 ]--;
    }

    // Copy the output array to arr[], so that arr[] now
    // contains sorted numbers according to current digit
    for (i = 0; i < n; i++)
        arr[i] = output[i];
}

// The main function to that sorts arr[] of size n using
// Radix Sort
template <typename T>
void radixsort(T* arr, int n)
{
    // Find the maximum number to know number of digits
    T m = getMax(arr, n);

    // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number
    for (unsigned int exp = 1; m/exp > 0; exp *= 10)
        countSort(arr, n, exp);
}


#endif //SFCSORTBENCH_RADIX_H
