#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <cstdio>
#include <cmath>
#include <regex>
#include <chrono>

#define threads 4

//CUDA
#include <cuda_runtime.h>
#include <cuda_profiler_api.h>

using namespace std;
using namespace std::chrono;

//TRANSKRYPCJA (DNA -> RNA)
__device__ void transcription(char* dna, char* rna, int index, int length_of_parts, int additional){  
  for(int i = 0; i < length_of_parts + additional; i++){
    //A-U
    if(dna[index * length_of_parts + i] == 'A'){
      rna[index * length_of_parts + i] = 'U';
    }
    //T-A
    else if(dna[index * length_of_parts + i] == 'T'){
      rna[index * length_of_parts + i] = 'A';
    } 
    //G-C
    else if(dna[index * length_of_parts + i] == 'G'){
      rna[index * length_of_parts + i] = 'C';
    }
    //C-G
    else if(dna[index * length_of_parts + i] == 'C'){
      rna[index * length_of_parts + i] = 'G';
    }
  }
}

__global__ void transform(char* dna, char* rna, char* result, int length_of_parts, int length, char* chain, char* code){
  int index = threadIdx.x;
  int additional = 0;

  if(index == threads - 1){
    additional = length - (length_of_parts * (threads - 1)) - length_of_parts;
  }
  transcription(dna, rna, index, length_of_parts, additional);
}

int main()
{
  string line, info;
  ifstream DNAFile("100.fasta");
  
  //ignorowanie pierwszej linijki, ponieważ są w niej zapisane dane genomu, a nie sam genom
  getline (DNAFile, info);

  //wczytywanie z pliku 
  string dna0, dna;
  while (getline (DNAFile, line)) {
    dna0 = dna0 + line;
  }
  for(char letter: dna0){
      dna.insert(0, string(1, letter));
  }

  char* charredDNA = (char*)malloc(sizeof(char) * dna.length());
  for(int i = 0; i < dna.length(); ++i){
    charredDNA[i] = dna[i];
  }
  charredDNA[dna.length()] = '\0';

  int length_of_parts = dna.length()/threads;

  char* cudaDNA;
  cudaMalloc((void**)&cudaDNA, sizeof(char) * dna.length());
  char* cudaResult;
  cudaMalloc((void**)&cudaResult, sizeof(char) * length_of_parts * threads);
  char* cudaRNA;
  cudaMalloc((void**)&cudaRNA, sizeof(char) * length_of_parts * threads);
  char* cudaChain;
  cudaMalloc((void**)&cudaChain, sizeof(char) * length_of_parts * threads / 3);
  char* cudaCode;
  cudaMalloc((void**)&cudaCode, sizeof(char) * threads * 3);
  cudaMemcpy(cudaDNA, charredDNA, sizeof(char) * dna.length(), cudaMemcpyHostToDevice);

  auto start = high_resolution_clock::now();
  transform<<<1, threads>>>(cudaDNA, cudaRNA, cudaResult, length_of_parts, dna.length(), cudaChain, cudaCode);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  //cout << "Transcription took " << duration.count() << " microseconds." << '\n';

  cudaFree(cudaDNA);
  free(charredDNA);

}
