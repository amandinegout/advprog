// ECE4122/6122 RSA Encryption/Decryption assignment
// Fall Semester 2015

#include <iostream>
#include "RSA_Algorithm.h"

using namespace std;

int main()
{
  // Instantiate the one and only RSA_Algorithm object
  RSA_Algorithm RSA;
  
  // Loop from sz = 32 to 1024 inclusive
  // for each size choose 10 different key pairs
  // For each key pair choose 10 differnt plaintext 
  // messages making sure it is smaller than n.
  // If not smaller then n then choose another
  // For eacm message encrypt it using the public key (n,e).
  // After encryption, decrypt the ciphertext using the private
  // key (n,d) and verify it matches the original message.


  for (size_t sz=32; sz<=1024; sz=sz*2){
	  //cout << "Size of p,q is " << sz << endl;
	  for (size_t i=0;i<10;i++){
		  RSA.GenerateRandomKeyPair(sz);
		  RSA.PrintNDE();
		  for (size_t j=0;j<10;j++){
			  mpz_class m = RSA.rng.get_z_bits(2*sz-1);
			  while (m>=RSA.n){
			 	  m = RSA.rng.get_z_bits(2*sz);
			  }
			  RSA.PrintM(m);
			  mpz_class c = RSA.Encrypt(m);
			  RSA.PrintC(c);
			  mpz_class m2 = RSA.Decrypt(c);
			  //RSA.PrintM(m2);
		  }
	  }
  }
  
return 0;
}
  
