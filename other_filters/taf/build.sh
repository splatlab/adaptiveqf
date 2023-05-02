g++ -std=c++11 src/taf.cc -o main src/arcd.o src/murmur3.o src/bit_util.o src/set.o -Wall -m64 -I. -Iinclude -Ofast -DNDEBUG -lpthread -lssl -lcrypto
