#include "LowMC.h"
#include <iostream>


//////////////////
//     MAIN     //
//////////////////

int main () {
    // Example usage of the LowMC class
    // Instantiate a LowMC cipher instance called cipher using the key '1'.
    LowMC cipher(1);
    block m = 0xFFD5;
    
    std::cout << "Plaintext:" << std::endl;
    std::cout << m << std::endl;
    m = cipher.encrypt( m );
    std::cout << "Ciphertext:" << std::endl;
    std::cout << m << std::endl;
    m = cipher.decrypt( m );
    std::cout << "Encryption followed by decryption of plaintext:" << std::endl;
    std::cout << m << std::endl;
   
    return 0;
}
