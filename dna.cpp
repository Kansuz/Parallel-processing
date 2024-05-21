#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <cstdio>
#include <cmath>

//OPENMP
#include <chrono> 
#include <omp.h>

using namespace std;

//DNA(3'->5') => RNA(5'->3')
string transcription(string dna){
  string rna;
  int len = dna.length();
  vector<string> acids(len);
  
  #pragma omp parallel for
  for(int i = 0; i < len; i++){

    //A-U
    if(dna[i] == 'A'){
        acids[i] = 'U';
    }
    //T-A
    else if(dna[i] == 'T'){
        acids[i] = 'A';
    } 
    //G-C
    else if(dna[i] == 'G'){
        acids[i] = 'C';
    }
    //C-G
    else if(dna[i] == 'C'){
        acids[i] = 'G';
    }
  }
  for (int v = 0; v < len; v++){
  rna = rna + acids[v];
  //printf("%d, %s\n", v, acids[v].c_str());
  }
  return rna;
}

//RNA(5'->3') => łańcuch białkowy
string translation(string rna){
  //UUU, UUC, UUA -> fenyloalanina
  //UUG, CUU, CUC, CUA, CUG -> leucyna
  //AUU, AUC, AUA -> izoleucyna
  //AUG -> START metionina
  //GUU, GUC, GUA, GUG -> walina
  //UCU, UCC, UCA, UCG, AGU, AGC -> seryna
  //CCU, CCC, CCA, CCG -> prolina
  //ACU, ACC, ACA, ACG -> treonina
  //GCU, GCC, GCA, GCG -> alanina
  //UAU, UAC -> tyrozyna
  //UAA, UAG, UGA -> STOP
  //CAU, CAC -> histydyna
  //CAA, CAG -> glutamina
  //AAU, AAC -> asparagina
  //AAA, AAG -> lizyna
  //GAU, GAC -> kwas asparaginowy
  //GAA, GAG -> kwas glutaminowy
  //UGU, UGC -> cysteina
  //UGG -> tryptofan 
  //CGU, CGC, CGA, CGG, AGA, AGG -> arginina
  //GGU, GGC, GGA, GGG -> glicyna
  string chain;
  if (rna.length()%3 == 1){
    rna.pop_back();
  }
  else if (rna.length()%3 == 2){
    rna.pop_back();
    rna.pop_back();
  }

  int len = rna.length()/3;
  vector<string> aminokwas(len);
  
  #pragma omp parallel for //num_threads(100)
  for(int i = 0; i < rna.length(); i=i+3){

    //int t = omp_get_thread_num();
    //printf("Hello from thread %d\n", t);
    
    string code = rna.substr(i, 3);
    string amino;
    //fenyloalanina
    if(code == "UUU" || code == "UUC" || code == "UUA"){
      amino =  "fenyloalanina-";
    }
    //leucyna
    else if(code == "UUG" || code == "CUU" || code == "CUC" || code == "CUA" || code == "CUG"){
      amino =  "leucyna-";
    }
    //izoleucyna
    else if(code == "AUU" || code == "AUC" || code == "AUA"){
      amino =  "izoleucyna-";
    }
    //metionina
    else if(code == "AUG"){
      amino =  "(START)metionina-";
    }
    //walina
    else if(code == "GUU" || code == "GUC" || code == "GUA" || code == "GUG"){
      amino =  "walina-";
    }
    //seryna
    else if(code == "UCU" || code == "UCC" || code == "UCA" || code == "UCG" || code == "AGU" || code == "AGC"){
      amino =  "seryna-";
    }
    //prolina
    else if(code == "CCU" || code == "CCC" || code == "CCA" || code == "CCG"){
      amino =  "prolina-";
    }
    //treonina
    else if(code == "ACU" || code == "ACC" || code == "ACA" || code == "ACG"){
      amino =  "treonina-";
    }
    //alanina
    else if(code == "GCU" || code == "GCC" || code == "GCA" || code == "GCG"){
      amino =  "alanina-";
    }
    //tyrozyna
    else if(code == "UAU" || code == "UAC"){
      amino =  "tyrozyna-";
    }
    //STOP
    else if(code == "UAA" || code == "UAG" || code == "UGA"){
      amino = "STOP";
    }
    //histydyna
    else if(code == "CAU" || code == "CAC"){
      amino =  "histydyna-";
    }
    //glutamina
    else if(code == "CAA" || code == "CAG"){
      amino =  "glutamina-";
    }
    //asparagina
    else if(code == "AAU" || code == "AAC"){
      amino =  "asparagina-";
    }
    //lizyna
    else if(code == "AAA" || code == "AAG"){
      amino =  "lizyna-";
    }
    //kwas asparginowy
    else if(code == "GAU" || code == "GAC"){
        amino =  "kwas asparginowy-";
    }
    //kwas glutaminowy
    else if(code == "GAA" || code == "GAG"){
      amino =  "kwas glutaminowy-";
    }
    //cysteina
    else if(code == "UGU" || code == "UGC"){
      amino =  "cysteina-";
    }
    //tryptofan
    else if(code == "UGG"){
      amino = "tryptofan-";
    }
    //arginina
    else if(code == "CGU" || code == "CGC" || code == "CGA" || code == "CGG" || code == "AGA" || code == "AGG"){
      amino = "arginina-";
    }
    //glicyna
    else if(code == "GGU" || code == "GGC" || code == "GGA" || code == "GGG"){
      amino = "glicyna-";
      //aminokwas[i] = "glicyna";
    }
    else{
      amino = "x";
    }
    aminokwas[i/3] = amino;
    //printf("%d, %s\n", i, amino.c_str());
  }
  for (int v = 0; v < len; v++){
    chain = chain + aminokwas[v];
  }

  return chain;
}

int main()
{
  string line;
  string trash;
  ifstream DNAFile("file.txt");
  
  //ignorowanie pierwszej linijki, ponieważ są w niej zapisane dane genomu, a nie sam genom
  getline (DNAFile, trash);
  // for(int i = 0; i < trash.length(); i++){
  //    cout << trash[i];
  // }

  //wczytywanie z pliku 
  string dna0, dna;
  while (getline (DNAFile, line)) {
    dna0 = dna0 + line; //jeżeli będzie odwrócone
    //dna = line + dna; //jeżeli nie będzie odwrócone
  }
  //jeżeli będzie odwrócone
  for(char letter: dna0){
      dna.insert(0, string(1, letter));
  }
  
  // UAC GGG GGG UAG
  //tyrozyna - glicyna - stop

  //PO POPRAWCE:
  //DNA: CTACCCCCCGTA
  //RNA: GAU GGG GGG CAU
  //A:C T:A C:C G:T

  //cout << "DNA:" << dna << '\n';

  string rna;
  cout << "DNA:" << dna;
  rna = transcription(dna);
  cout << "RNA:" << rna << '\n';

  string chain;
  chain = translation(rna);
  cout << chain << '\n';

  ofstream result("result.txt");
  result << chain;
  result.close();

  //statystyki--------------------------
  //1. Jaka zasada najczęściej występuje po danej zasadzie
  //A T C G
  vector<int> A = {0, 0, 0, 0};
  vector<int> T = {0, 0, 0, 0};
  vector<int> C = {0, 0, 0, 0};
  vector<int> G = {0, 0, 0, 0};
  for (int i = 0; i < dna.length(); i++){
    if (dna[i] == 'A'){
      if (dna[i+1] == 'A'){
        A[0] += 1;
      }
      else if (dna[i+1] == 'T'){
        A[1] += 1;
      }
      else if (dna[i+1] == 'C'){
        A[2] += 1;
      }
      else if (dna[i+1] == 'G'){
        A[3] += 1;
      }
    }
    else if (dna[i] == 'T'){
      if (dna[i+1] == 'A'){
        T[0] += 1;
      }
      else if (dna[i+1] == 'T'){
        T[1] += 1;
      }
      else if (dna[i+1] == 'C'){
        T[2] += 1;
      }
      else if (dna[i+1] == 'G'){
        T[3] += 1;
      }
    }
    else if (dna[i] == 'C'){
      if (dna[i+1] == 'A'){
        C[0] += 1;
      }
      else if (dna[i+1] == 'T'){
        C[1] += 1;
      }
      else if (dna[i+1] == 'C'){
        C[2] += 1;
      }
      else if (dna[i+1] == 'G'){
        C[3] += 1;
      }
    }
    else if (dna[i] == 'G'){
      if (dna[i+1] == 'A'){
        G[0] += 1;
      }
      else if (dna[i+1] == 'T'){
        G[1] += 1;
      }
      else if (dna[i+1] == 'C'){
        G[2] += 1;
      }
      else if (dna[i+1] == 'G'){
        G[3] += 1;
      }
    }
  }
  cout << "DNA:" << dna << '\n';
  for (int i = 0; i < 4; i++){
    printf("%d", A[i]);
  }
  cout << "\n";
  for (int i = 0; i < 4; i++){
    printf("%d", T[i]);
  }
  cout << "\n";
  for (int i = 0; i < 4; i++){
    printf("%d", C[i]);
  }
  cout << "\n";
  for (int i = 0; i < 4; i++){
    printf("%d", G[i]);
  }
  cout << "\n";

  //DNA: CTACCCCCCGTA
  //A 0010 -> C
  //T 2000 -> A
  //C 0150 -> C
  //G 0100 -> T
     
}

