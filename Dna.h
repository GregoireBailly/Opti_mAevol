//
// Created by arrouan on 01/10/18.
//

#pragma once

#include <cstdio>
#include <cstring>
#include <cassert>
#include <cstdint>
#include <vector>
#include <zlib.h>

#include "Threefry.h"

static unsigned char or_mask[] = {0b00000001, 0b00000010, 0b00000100, 0b00001000, 
                  0b00010000, 0b00100000, 0b01000000, 0b10000000};
static unsigned char and_mask[] ={0b11111110, 0b11111101, 0b11111011, 0b11110111, 
                  0b11101111, 0b11011111, 0b10111111, 0b01111111};

class Bitset {
    public:
        Bitset() = default;

        Bitset(const Bitset& clone){
            data_ = new unsigned char [clone.length_];
            length_ = clone.length_;
            max_length_ = length_;
            strcpy(data_, clone.data_);
        }

        Bitset(char* genome, int length){
          data_ = new unsigned char [length];
          length_ = length;
          max_length_ = length_;
          strcpy(data_, genome);
        }

        Bitset(int length) {
          data_ = new unsigned char [length];
          length_ = length;
          max_length_ = length_;
        }

        ~Bitset() = default;

        bool get(int pos){
            return ((data_[pos/8] & or_mask[pos%8]) != 0);
        }

        void flip(int pos){
            data_[pos/8] = data_[pos/8] ^ or_mask[pos%8];
        }

        void set(int pos){
            data_[pos/8] = data_[pos/8] | or_mask[pos%8];
        }

        void clear(int pos){
            data_[pos/8] = data_[pos/8] & and_mask[pos%8];
        }

        void remove(int pos1, int pos2){
            for(int i=1; i<length_-pos2; i++){
                if(get(pos2+i))
                    set(pos1+i);
                else
                    clear(pos1+i);
            }
            for(int i=pos1+length_-pos2; i<length_; i++){
                clear(i);
            }
        }

        void insert(Bitset add, int pos){
          if(length_+add.length_ > max_length_) {
            max_length_ = max_length_ * 2;
            unsigned char* data_tmp = new unsigned char [max_length_];
            strcpy(data_tmp, data_);
            data_ = data_tmp;
          }
          for(int i=length_-1; i>=pos; i--)
          {
              if(get(i))
                  set(add.length_+i);
              else
                  clear(add.length_+i);
          }
          for(int i=0; i<add.length_; i++)
          {
              if(add.get(i))
                  set(pos+i);
              else
                  clear(pos+i);
          }
        }

    private:
        unsigned char* data_;
        int length_;
        int max_length_;
};

constexpr int8_t CODON_SIZE = 3;

constexpr const char* PROM_SEQ = "0101011001110010010110";
constexpr const char* SHINE_DAL_SEQ = "011011000";
constexpr const char* PROTEIN_END = "001"; // CODON_STOP

class ExpManager;

class Dna {

 public:
  Dna() = default;

  Dna(const Dna& clone);

  Dna(int length, Threefry::Gen& rng);

  Dna(char* genome, int length);

  Dna(int length);

  ~Dna() = default;

  int length() const;

  void save(gzFile backup_file);
  void load(gzFile backup_file);

  void set(int pos, char c);

  /// Remove the DNA inbetween pos_1 and pos_2
  void remove(int pos_1, int pos_2);

  /// Insert a sequence of a given length at a given position into the DNA of the Organism
  void insert(int pos, Bitset & seq);

  /// Insert a sequence of a given length at a given position into the DNA of the Organism
  void insert(int pos, Dna* seq);

  void do_switch(int pos);

  void do_duplication(int pos_1, int pos_2, int pos_3);

  int promoter_at(int pos);

  int terminator_at(int pos);

  bool shine_dal_start(int pos);

  bool protein_stop(int pos);

  int codon_at(int pos);

  Bitset seq_;
  int32_t length_;
};
