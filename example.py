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
