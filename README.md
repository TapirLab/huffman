# Huffman Coding Python & MATLAB Implementations

![Tapir Lab.](http://tapirlab.com/wp-content/uploads/2020/10/tapir_logo.png)
This repository consists of `MATLAB` and `Python` implementations of Huffman coding. Huffman source coding is introduced by **David Albert Huffman** and published with the title of "[A Method for the Construction of Minimum-Redundancy Codes](https://ieeexplore.ieee.org/document/4051119)" in Proceedings of I.R.E., September 1952. 

## Description
Huffman coding is a minimum-redundancy and variable-length source coding method. These terms require explanation before proceeding further. First of all, minimum-redundancy means that coding is performed in such a way that average bit length is minimized. Minimizing the average bit length means that maximum compression (optimum coding) is achieved without any loss of information. Optimum coding is performed by assigning the shortest code to the most probable symbol and the longest code to the least probable symbol. This procedure introduces the variable-length coding concept. Additionally, Huffman is a prefix code. Code initials are the distinguishing property of each code. That means, each digit of code represents either left or right in a binary Huffman tree, and code words are located only on the leaves. For instance, assume that you have the text `AAAAAAAAAAAABBBBBBCCCDDD` and want to encode it. Observe that text consists of four unique symbols `['A', 'B', 'C', 'D']` and their frequencies in the text is `[0.50, 0.25, 0.125, 0.125]` respectively. If you do not consider the frequency of symbols, each can be coded by using two bits (`A`->`00`, `B`->`01`, `C`->`11`, `D`->`10`). However, additional information about the frequency of symbols makes more efficient coding possible. By assigning `0` to `A`, `10` to `B`, `110` to `C`, and `111` to D, Huffman coding can be achieved. Initially, using three bits for coding `C` and `D` may seem ineffective but if you multiply the frequency of each symbol with code length and sum them up (averaging) you will get `1.75` which is equals to Shannon's Entropy. In order to apply Huffman coding, the input should have a finite number of symbols and produced Huffman codes restricted by the followings:
1. Each code is assigned to one and only one symbol.
2. A code can not start with any other codes i.e. each code is unique.
3. For binary coding, there can be only two different symbols that have the same bit sequences except for the last digit.
4. Let N represent the length of the code. Each possible N-1 bit sequence must be used as a message code or prefix to another message code.

## Prerequisites

* Python3
* numpy

**Note**: `1.19.4` version of `numpy` has some issues on Windows. Thus, `requirements.txt` specifies version as `1.19.3`. Please pay attention to this detail if you want to use on Windows.

## Folder Structure

```
huffman
|── matlab_implementation
    |── Source codes of Huffman coding Matlab implementation
|── sample_input
    |── An example text file
|── sample_output
    |── If reports are generated, they will be located under this folder
|── example.py
|── huffman.py
|── LICENSE
|── README.md
|── requirements.txt
```
## Example 

After installing requirement (numpy), this repo can be tested by executing the `example.py` script. It will create a random text which includes 10 unique upper case ASCII characters and has a length of 500. Then, a code dictionary for this input text will be created. After creation of dictionary, text will be encoded and decoded. Reports will be saved under `sample_output` folder. Content of `example.py` can be seen below.

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file shows how Huffman coding and decoding are performed. Initially,
a random text is generated. Then, Huffman coding is performed, and results are
saved as different files to the same folder. For more information please read
README.md or visit GitHub page

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Huffman Coding
%% -------------------
%% $Author: Halil Said Cankurtaran$,
%% $Date: August 8th, 2022$,
%% $Revision: 1.2$
%% $Tapir Lab.$
%% $Copyright: Tapir Lab.$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
import os
import huffman as huff

if __name__ == "__main__":
    # If you want to read input text directly from a file please configure followings
    # and comment out random text generation
    # file_name = 'input_text.txt'
    # path_to_text_file = os.path.join('.', 'sample_input', file_name)
    # input_text = huff.read_text_from_file(path_to_text_file)

    # Initilize required variables to create random text
    type_of_distribution = 'exponential'
    number_of_realization = 500
    number_of_symbols = 10
    decay_rate = 0.01
    # Generates different random upper case ASCII characters as specified by `number_of_symbols`
    random_symbols = huff.generate_n_random_symbols(number_of_symbols)
    freqs = huff.generate_random_symbol_frequencies(type_of_distribution,
                                               number_of_realization,
                                               number_of_symbols,
                                               decay_rate)
    input_text = huff.generate_random_text_with_predefined_frequencies(freqs, random_symbols, 200)

    # Calculate frequency of symbols
    symbols_and_freqs, unique_symbols = huff.calculate_frequency_of_symbols_from_text(input_text)
    # Build Tree
    tree = huff.build_tree(symbols_and_freqs)
    # Traverse tree and build dictionary
    # `tree_type=0` means, left -> 0 and right -> 1. `tree_type=1` means, left -> 1 and right -> 0.
    dictionary = huff.traverse_tree(tree, unique_symbols, tree_type=0)

    # Encode the input text based on created dictionary
    encoded_text = huff.encode_message(dictionary, input_text)
    # Decode the encoded text based on created dictionary 
    decoded_text = huff.decode_message(encoded_text, dictionary)
    
    # Create reports
    huff.produce_huffman_report(symbols_and_freqs, dictionary)
    huff.huffman_coded_string_report(input_text, encoded_text, symbols_and_freqs, dictionary)
    huff.huffman_decoded_string_report(encoded_text, decoded_text, symbols_and_freqs, dictionary)
    
    # Calculation of average symbol length and Shannon entropy
    expected_length, shannon_entropy = huff.calculate_expected_lengh_and_shannon_entropy(symbols_and_freqs, dictionary)
    # Print Results
    print("A random text has been created and encoded. Dictionary can be seen below.")
    print("--------------------------")
    print("Random text: ", input_text)
    print("--------------------------")
    print("Symbols and their frequencies: ", symbols_and_freqs)
    print("--------------------------")
    print("Dictionary: ", dictionary)
    print("--------------------------")
    print("Encoded Text: ", encoded_text)
    print("--------------------------")
    print("Decoded Text: ", decoded_text)
    print("--------------------------")
    print(f'Expected length = {expected_length}')
    print(f'Shannon Entropy = {shannon_entropy}')
    print("--------------------------")
    print("Reports also saved into sample_output folder")
```

## License

The software is licensed under the MIT License.
