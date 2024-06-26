#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <cstdio>
#include <cmath>
#include <regex>
#include <chrono>

using namespace std;
using namespace std::chrono;

//TRANSKRYPCJA (DNA -> RNA)
string transcription(string dna){
  string rna;
  int len = dna.length();
  vector<string> acids(len);
  
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
  }
  return rna;
}

//SPLICING
string splicing(string rna){
  string spliced_rna;
  regex intron("GU.*?AG");
  if (regex_search(rna, intron)){
    spliced_rna = regex_replace(rna, intron, "");
  }
  return spliced_rna;
}

//TRANSLACJA (RNA -> łańcuch białkowy)
string translation(string rna){
  string chain;

  int len = rna.length();
  vector<string> aminokwas(len);
  
  for(int i = 0; i < len; i=i+3){

    string code = rna.substr(i, 3);

    //fenyloalanina
    if(code == "UUU" || code == "UUC" || code == "UUA"){
      aminokwas[i] =  "fenyloalanina-";
    }
    //leucyna
    else if(code == "UUG" || code == "CUU" || code == "CUC" || code == "CUA" || code == "CUG"){
      aminokwas[i] =  "leucyna-";
    }
    //izoleucyna
    else if(code == "AUU" || code == "AUC" || code == "AUA"){
      aminokwas[i] =  "izoleucyna-";
    }
    //metionina
    else if(code == "AUG"){
      aminokwas[i] =  "(START)metionina-";
    }
    //walina
    else if(code == "GUU" || code == "GUC" || code == "GUA" || code == "GUG"){
      aminokwas[i] =  "walina-";
    }
    //seryna
    else if(code == "UCU" || code == "UCC" || code == "UCA" || code == "UCG" || code == "AGU" || code == "AGC"){
      aminokwas[i] =  "seryna-";
    }
    //prolina
    else if(code == "CCU" || code == "CCC" || code == "CCA" || code == "CCG"){
      aminokwas[i] =  "prolina-";
    }
    //treonina
    else if(code == "ACU" || code == "ACC" || code == "ACA" || code == "ACG"){
      aminokwas[i] =  "treonina-";
    }
    //alanina
    else if(code == "GCU" || code == "GCC" || code == "GCA" || code == "GCG"){
      aminokwas[i] =  "alanina-";
    }
    //tyrozyna
    else if(code == "UAU" || code == "UAC"){
      aminokwas[i] =  "tyrozyna-";
    }
    //STOP
    else if(code == "UAA" || code == "UAG" || code == "UGA"){
      aminokwas[i] = "STOP";
    }
    //histydyna
    else if(code == "CAU" || code == "CAC"){
      aminokwas[i] =  "histydyna-";
    }
    //glutamina
    else if(code == "CAA" || code == "CAG"){
      aminokwas[i] =  "glutamina-";
    }
    //asparagina
    else if(code == "AAU" || code == "AAC"){
      aminokwas[i] =  "asparagina-";
    }
    //lizyna
    else if(code == "AAA" || code == "AAG"){
      aminokwas[i] =  "lizyna-";
    }
    //kwas asparginowy
    else if(code == "GAU" || code == "GAC"){
      aminokwas[i] =  "kwas asparginowy-";
    }
    //kwas glutaminowy
    else if(code == "GAA" || code == "GAG"){
      aminokwas[i] =  "kwas glutaminowy-";
    }
    //cysteina
    else if(code == "UGU" || code == "UGC"){
      aminokwas[i] =  "cysteina-";
    }
    //tryptofan
    else if(code == "UGG"){
      aminokwas[i] = "tryptofan-";
    }
    //arginina
    else if(code == "CGU" || code == "CGC" || code == "CGA" || code == "CGG" || code == "AGA" || code == "AGG"){
      aminokwas[i] = "arginina-";
    }
    //glicyna
    else if(code == "GGU" || code == "GGC" || code == "GGA" || code == "GGG"){
      aminokwas[i] = "glicyna-";
    }
    else{
      aminokwas[i] = "x";
    }
  }
  for (int v = 0; v < len; v++){
    chain = chain + aminokwas[v];
  }
  return chain;
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
  cout << "DNA: " << dna << '\n';

  int threads = 4;

  //transkrypcja
  string rna;
  auto start_rna = high_resolution_clock::now();
  rna = transcription(dna);
  auto stop_rna = high_resolution_clock::now();
  cout << "RNA: " << rna << '\n';
  auto duration_rna = duration_cast<microseconds>(stop_rna - start_rna);

  //splicing
  string spliced_rna;
  spliced_rna = splicing(rna);
  cout << "Spliced RNA: " << spliced_rna << '\n';

  //translacja
  string chain;
  auto start_chain = high_resolution_clock::now();
  chain = translation(spliced_rna);
  auto stop_chain = high_resolution_clock::now();
  cout << "Amino acid sequence: " << chain << '\n';
  auto duration_chain = duration_cast<microseconds>(stop_chain - start_chain);

  cout << "Time taken by transcription: " << duration_rna.count() << " microseconds" << '\n';
  cout << "Time taken by translation: " << duration_chain.count() << " microseconds" << '\n';
  
  ofstream result("result_seq.txt");
  result << chain;
  result.close();
  }