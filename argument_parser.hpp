/**
 * shark - Mapping-free filtering of useless RNA-Seq reads
 * Copyright (C) 2019 Tamara Ceccato, Luca Denti, Yuri Pirola, Marco Previtali
 *
 * This file is part of shark.
 *
 * shark is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * shark is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with shark; see the file LICENSE. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#ifndef _ARGUMENT_PARSER_HPP_
#define _ARGUMENT_PARSER_HPP_

#include <iostream>
#include <sstream>
#include <getopt.h>

static const char *USAGE_MESSAGE =
"Usage: shark -r <references> -1 <sample1> [OPTIONAL ARGUMENTS]\n"
"\n"
"Arguments:\n"
"      -r, --reference                   reference sequences in FASTA format (can be gzipped)\n"
"      -1, --sample1                     sample in FASTQ (can be gzipped)\n"
"\n"
"Optional arguments:\n"
"      -h, --help                        display this help and exit\n"
"      -2, --sample2                     second sample in FASTQ (optional, can be gzipped)\n"
"      -o, --out1                        first output sample in FASTQ (default: sharked_sample.1)\n"
"      -p, --out2                        second output sample in FASTQ (default: sharked_sample.2)\n"
"      -k, --kmer-size                   size of the kmers to index (default:17, max:31)\n"
"      -c, --confidence                  confidence for associating a read to a gene (default:0.6)\n"
"      -b, --bf-size                     bloom filter size in Kb (default:1024)\n"
"      -q, --min-base-quality            minimum base quality (assume FASTQ Illumina 1.8+ Phred scale, default:0, i.e., no filtering)\n"
"      -s, --single                      report an association only if a single gene is found\n"
"      -t, --threads                     number of threads (default:1)\n"
"      -m, --method                      subject of the condition [base / kmer] (default: base)\n"
"      -x, --xxhash                      number of hash functions\n"
"      -v, --verbose                     verbose mode\n";

namespace opt {
  static std::string fasta_path = "";
  static std::string sample1_path = "";
  static std::string sample2_path = "";
  static std::string out1_path = "";
  static std::string out2_path = "";
  static bool paired_flag = false;
  static uint k = 17;
  static double c = 0.6;
  static uint64_t bf_size = ((uint64_t)0b1 << 20);
  static char min_quality = 0;
  static bool single = false;
  static std::string method = "";
  static int nHash = 1;
  static bool verbose = false;
  static int nThreads = 1;
}

static const char *shortopts = "t:r:1:2:o:p:k:c:b:q:m:x:y:svh";

static const struct option longopts[] = {
  {"reference", required_argument, NULL, 'r'},
  {"threads", required_argument, NULL, 't'},
  {"sample1", required_argument, NULL, '1'},
  {"sample2", required_argument, NULL, '2'},
  {"out1", required_argument, NULL, 'o'},
  {"out2", required_argument, NULL, 'p'},
  {"kmer-size", required_argument, NULL, 'k'},
  {"confidence", required_argument, NULL, 'c'},
  {"bf-size", required_argument, NULL, 'b'},
  {"min-base-quality", required_argument, NULL, 'q'},
  {"single", no_argument, NULL, 's'},
  {"method", required_argument, NULL, 'm'},
  {"xxhash", required_argument, NULL, 'x'},
  {"verbose", no_argument, NULL, 'v'},
  {"help", no_argument, NULL, 'h'},
  {NULL, 0, NULL, 0}
};

void parse_arguments(int argc, char **argv) {
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1; ) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'r':
      arg >> opt::fasta_path;
      break;
    case 't':
      arg >> opt::nThreads;
      if(opt::nThreads <= 0) {
        std::cerr << "USAGE_MESSAGE";
        std::cerr << "shark: at least 1 thread is required." << std::endl
                  << "aborting..." << std::endl;
        exit(EXIT_FAILURE);
      }
      break;
    case '1':
      arg >> opt::sample1_path;
      break;
    case '2':
      arg >> opt::sample2_path;
      opt::paired_flag = true;
      break;
    case 'o':
      arg >> opt::out1_path;
      break;
    case 'p':
      arg >> opt::out2_path;
      break;
    case 'k':
      arg >> opt::k;
      if(opt::k == 0 or opt::k > 31) {
        std::cerr << USAGE_MESSAGE;
        std::cerr << "shark: k must be in the range [1, 31]." << std::endl
                  << "aborting..." << std::endl;
        exit(EXIT_FAILURE);
      }
      break;
    case 'c':
      arg >> opt::c;
      if(opt::c < 0 or opt::c > 1) {
        std::cerr << "shark: c must be in the range [0, 1]." << std::endl
                  << "aborting..." << std::endl;
        exit(EXIT_FAILURE);
      }
      break;
    case 'b':
      // Let's consider this as Kb
      arg >> opt::bf_size;
      opt::bf_size = opt::bf_size * ((uint64_t)0b1 << 10);
      break;
    case 'q':
      int mq;
      arg >> mq;
      if(mq < 0) {
        std::cerr << USAGE_MESSAGE;
        std::cerr << "shark: q must be a positive value." << std::endl
                  << "aborting..." << std::endl;
        exit(EXIT_FAILURE);
      }
      opt::min_quality = static_cast<char>(mq);
      break;
    case 's':
      opt::single = true;
      break;
    case 'm':
      arg >> opt::method;
      break;
    case 'x':
      arg >> opt::nHash;
      if(opt::nHash <= 0) {
        std::cerr << "USAGE_MESSAGE";
        std::cerr << "shark: at least 1 hash function is required." << std::endl
                  << "aborting..." << std::endl;
        exit(EXIT_FAILURE);
      }
	  break;
    case 'v':
      opt::verbose = true;
      break;
    case 'h':
      std::cerr << USAGE_MESSAGE;
      exit(EXIT_SUCCESS);
    default:
      std::cerr << "shark : unknown argument" << std::endl;
      std::cerr << "\n" << USAGE_MESSAGE;
      exit(EXIT_FAILURE);
    }
  }

  if (opt::fasta_path == "" || opt::sample1_path == "") {
    std::cerr << "shark : missing required files" << std::endl;
    std::cerr << "\n" << USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }

  if(opt::out1_path == "") {
    opt::out1_path = "sharked_sample.1";
  }
  if(opt::out2_path == "" and opt::sample2_path != "") {
    opt::out2_path = "sharked_sample.2";
  }
}

#endif
