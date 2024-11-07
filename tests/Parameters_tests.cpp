#include "../src/parameters.h"
#include <iostream>

int main() {
    // Create a Parameters object
    Parameters params;

    // Test storing different types of values
    params["doubleValue"] = Variable(3.14);
    params["intValue"] = Variable(42);
    params["boolValue"] = Variable(true);
    params["stringValue"] = Variable(std::string("Hello, World!"));
    params["vectorValue"] = Variable(std::vector<double>{1.1, 2.2, 3.3});

    // Retrieve and print each type to test the get() method and implicit conversions
    std::cout << "doubleValue: " << params["doubleValue"].get<double>() << std::endl;
    std::cout << "intValue: " << params["intValue"].get<int>() << std::endl;
    std::cout << "boolValue: " << params["boolValue"].get<bool>() << std::endl;
    std::cout << "stringValue: " << params["stringValue"].get<std::string>() << std::endl;

    // Test vector retrieval
    std::vector<double> vec = params["vectorValue"].get<std::vector<double>>();
    std::cout << "vectorValue: ";
    for (double v : vec) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    // Test error handling by attempting to get the wrong type
    try {
        std::cout << "Attempting to retrieve intValue as a double: ";
        double incorrectType = params["intValue"].get<double>();
        std::cout << incorrectType << std::endl;
    } catch (const std::bad_variant_access&) {
        std::cerr << "Caught std::bad_variant_access: type mismatch!" << std::endl;
    }

    return 0;
}
