#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python implementation of Huffman Coding. This file includes functions that
can be used to perform a Huffman coding and decoding.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Huffman Coding
%% -------------------
%% $Author: Halil Said Cankurtaran$,
%% $Date: January 5th, 2020$,
%% $Revision: 1.1$
%% $Tapir Lab.$
%% $Copyright: Tapir Lab.$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
import numpy as np
import os

def read_text_from_file(path_to_text_file):
    """Reads the text file and returns whole content as a string.

    Args:
        path_to_text_file (str)

    Returns:
        A string that concatenates all the lines in a single line and uses 
        '\n' as the line separator. That means if the text file includes any '\n'
        it is also Huffman Coded.

        str: content of text file in a single line
    """
    with open(path_to_text_file, "r") as file:
        array_of_lines = [line.strip('\n') for line in file.readlines()]

    line_seperator = "\n"  # To join all the read lines into one string
    text = line_seperator.join(array_of_lines)  # Join into one string

    return text


def generate_random_symbol_frequencies(distribution_type='exponential',
                                       number_of_realization=50,
                                       number_of_symbols=6,
                                       decay_rate=0.01):
    """Generates random symbol frequencies based on parameters.

    Function returns the random frequency of occurrences based on the parameters.
    It is possible to generate uniformly or exponentially distributed symbol
    sequence frequencies for now. The function returns an exponentially distributed
    freqs with 50 realizations, 0.01 decay rate, and 6 symbols by default. For a
    better approximation of distribution, increase the number of realization value.
    If 'uniform', the decay rate has no effect.

    Args:
        distribution_type (str, optional): ['uniform', 0] or ['exponential', 1]. Defaults to 'exponential'.
        number_of_realization (int, optional): `size` parameter of random function. Defaults to 50.
        number_of_symbols (int, optional): Defaults to 6.
        decay_rate (float, optional): Decay rate of exp. dist. Defaults to 0.01.

    Returns:
        List: A list of numbers whose length equals to `number_of_symbols` and adds up to 1.
    """

    # Create random variables based on specified distribution
    if distribution_type in ['uniform', 0]:
        rands = np.random.uniform(0, number_of_symbols, number_of_realization)
    if distribution_type in ['exponential', 1]:
        rands = np.random.exponential(decay_rate, number_of_realization)
    else:
        raise ValueError('Specify istribution_type as \'uniform\' or \'exponential\'')

    # Create number of occurance by drawing histogram whose number of bins
    # equals to the number_of_symbols
    number_of_occurances, bins = np.histogram(rands, number_of_symbols)
    freq_list = number_of_occurances/np.sum(number_of_occurances)

    return sorted(freq_list, reverse=True)

def generate_n_random_symbols(n):
    """Returns `n` random upper case ASCII symbols

    Args:
        number_of_symbols (int): Number of symbols required

    Raises:
        ValueError: If n is greater than 26 and less than 0
    """
    if n in range(27):
        upper_case_ascii_characters = [chr(i) for i in range(65,91)]
        random_symbols = np.random.choice(upper_case_ascii_characters, n, replace=False)
    else:
        raise ValueError('n should be greater than zero and less than 26!')

    return random_symbols

def generate_random_text_with_predefined_frequencies(frequency_list,
                                                     symbols_list,
                                                     number_of_realization,
                                                     scramble=True):
    """Generate a random text with specified frequency list and symbols

    Args:
        frequency_list (list/np.array): List of frequencies corresponding to symbols.
        symbols_list ([type]): List of symbols.
        number_of_realization ([type]): Length of random text
        scramble (bool, optional): If true shuffles generated random text. Defaults to True.

    Returns:
        string: A randomly generated text based on given inputs.
    """

    if len(frequency_list) != len(symbols_list):
        raise ValueError('For each symbol, there must be a corresponding frequency. Lengths of frequency_list and symbol_list are different')

    tmp = np.array([]) # A temporary np array
    # Match symbols ith frequencies and iterate over
    for (symbol, freq) in zip(symbols_list, frequency_list):
        repetition = np.repeat(symbol,round(freq*number_of_realization))
        tmp = np.hstack((tmp, repetition))

    # If scramble is true, shuffle the tmp array
    if scramble:
        np.random.shuffle(tmp)
    
    # join list of symbols to string
    random_text = ""
    for elem in tmp:
        random_text += "".join(elem)

    return random_text


def calculate_frequency_of_symbols_from_text(text):
    """Finds unique symbols in the text and calculates their frequency

    Args:
        text (str): a string with symbols

    Returns:
        list: frequency of corresponding symbols
        set: unique symbols in the text
    """
    unique_symbols_in_text = set(text)  # Find each unique symbol in the text
    length_of_text = len(text)  # Store length of text to calculate frequency
    symbols_and_freqs = []  # Store probability of each symbol in text
    for symbol in unique_symbols_in_text:
        symbols_and_freqs.append([symbol, text.count(symbol)/length_of_text])

    # Sort frequencies in decreasing order
    symbols_and_freqs = sorted(symbols_and_freqs, key=lambda elem: elem[1], reverse=True)

    return symbols_and_freqs, unique_symbols_in_text

def check_order_of_symbols_and_freqs_list(symbols_and_freqs_list):
    """Checks whether provided symbols_freq_list is in decreasing order or not."""
    assume_max_freq = symbols_and_freqs_list[0][1]
    for elem in symbols_and_freqs_list:
        if elem[1] > assume_max_freq:
            raise ValueError('Provided list is not in decreasing order!')
    
    return True

def build_tree(symbols_and_freqs_in_decreasing_order):
    """Creates Huffman code tree of given frequency list which is ordered decreasingly.

    Steps in Huffman tree building:
        1. Pop least two elements of the tree,
        2. Sum them up and insert the result into the list again,
        3. Sort the list in decreasing order again
        4. Go to 1 while there are at least two elements in the list

    Args:
        symbols_and_freqs_in_decreasing_order (list): symbols and corresponding freqs
            As the loop proceeds, symbols are concatenated and freqs are summed

    Returns:
        dict: each node with child nodes such as: (freq_of_a > freq_of_b)
            {'ab': [sum_of_freqs, ['a', freq_of_a], ['b', freq_of_b]]}
    """
    # Check whether list is actually in decreasing order or not
    # Do not continue process in case of failure
    check_order_of_symbols_and_freqs_list(symbols_and_freqs_in_decreasing_order)
    tree = {}
    tmp = symbols_and_freqs_in_decreasing_order  # assigned to shorter variable
    while len(tmp) > 1:
        a, b = tmp[-2:]  # Get least two elements of array
        tmp = tmp[:-2]  # Remove least two elements of array
        string = a[0] + b[0]  # Store symbols of least two
        sum_of_symbols = a[1] + b[1]  # Sum freqs of least two elements
        tmp.append([string, sum_of_symbols])  # Append new node to the end
        # Sort array in decreasing order again
        tmp = sorted(tmp, key=lambda elem: elem[1], reverse=True)
        node = [sum_of_symbols, a, b]  # Create newly calculated node
        tree[string] = node  # Add node to the tree

    return tree


def build_dictionary(tree, dictionary, current_node, code, tree_type=0):
    """Traverses the tree and builds the code dictionary. Must be used in 
    `traverse_tree` function.

    The tree must be built by the `build_tree` function or in the same form.
    Performs depth-first search operation recursively, builds code dictionary.
    Assigns `0` to the left node and `1` to the right node by default. If you want to
    change this to the exact opposite, use `tree_type=1`

    Args:
        tree (dict): Code tree that is built with `build_tree` function
        dictionary (dict): Code dictionary, initially an empty dict
        current_node (str): The node which the function currently tests
        code (list): 

    Returns:
        dict: Huffman code dictionary
    """
    if tree_type == 1:
        codes = range(1, -1, -1)
    else: 
        codes = range(2)

    if current_node in list(tree.keys()):
        for (child, i) in zip(tree[current_node][1:], codes):
            code.append(f'{i}')
            build_dictionary(tree, dictionary, child[0], code, tree_type)
            code.pop()
    else:
        code_str = "".join(code)
        dictionary[current_node] = code_str

    return dictionary


def traverse_tree(tree, unique_symbols_in_text, tree_type=0):
    """Traverses created the tree and builds Huffman code dictionary

    A depth-first search is performed to the constructed tree. While traversing
    Assigns `0` to the left node and `1` to the right node by default. If you want to
    change this to the exact opposite, use `tree_type=1`

    Args:
        tree (dict): Code tree that is built with `build_tree` function
        unique_symbols_in_text (set): A set of unique symbols in text
        tree_type (int): 0 to assign ([Left, Right]) [0,1], 1 to assign [1,0]

    Returns:
        dict: Huffman code dictionary
    """
    code = []
    dictionary = {}
    tree_keys = list(tree.keys())
    current_node = tree_keys[-1]
    code_dictionary = build_dictionary(tree, dictionary, current_node, code, tree_type)

    return code_dictionary


def encode_message(dictionary, text):
    """Encodes given text according to the dictionary
    
    If symbols in the text are not a subset of symbols in the dictionary then,
    the function raises an error. Using a large dictionary instead of text-specific
    ones may result in inefficient coding!

    Args:
        dictionary (dict): Huffman code dictionary,
        test (str): text to encode according to a given dictionary

    Returns:
        str: encoded text

    Raises:
        Inconsistency error, if symbols in the text is not a subset of symbols in the dictionary
    """
    # In case a symbol is detected out of dictionary stop encoding!
    # Explain procedure 
    assert set(text).issubset(set(dictionary.keys())), 'A symbol detected out of dictionary!'
    for symbol in dictionary.keys():
        text = text.replace(symbol, dictionary[symbol])
    return text


def decode_message(encoded_text, dictionary):
    """Decodes encoded text according to the given dictionary.

    Huffman encoded messages can be decoded only if we have the original dictionary
    that is used to encode the message. If the encoded message includes an unknown code
    function raises an error.

    Args:
        encoded_text (str): Huffman encoded text
        dictionary (dict): A dict that includes corresponding symbols and codes
                           For example: {'a': '[code_a]', 'b': '[code_b]'}

    Returns:
        str: Decoded text according to given dictionary
    """
    decoded_text = ''
    # Since dictionary key is the "symbol" and value is the "code"
    # Keys and values must be swapped to access symbol based on read code
    reverse_dict = dict([(value,key) for (key,value) in dictionary.items()])
    longest_code = max([len(code) for code in reverse_dict.keys()])
    counter = 0
    while len(encoded_text) != 0:
        found = False
        while not found:
            tmp = encoded_text[:counter]
            assert (len(tmp) <= longest_code), f'{tmp} is not in the dictionary! Can not decode the message!'
            if tmp in reverse_dict.keys():
                decoded_text += reverse_dict[tmp]
                encoded_text = encoded_text[counter:]
                found = True
                counter = 0
            counter +=1

    return decoded_text

def calculate_expected_lengh_and_shannon_entropy(symbols_and_freqs, dictionary):
    """Calculate expected length of resulting coding process and Shannon entropy"""
    # merge codes with symbols and frequency lise
    symbols_and_freqs_codes = [elem + [dictionary[elem[0]]] for elem in symbols_and_freqs]
    # Extract freqs from smybols_freqs
    freqs = [elem[1] for elem in symbols_and_freqs]
    # Calculation of sum of multiple of frequency and code length - Expected Length
    expected_length = np.sum([freq*len(code) for (symbol, freq, code) in symbols_and_freqs_codes])
    # Calculation of sum of multiple of frequency and log2(frequency) - Shannon Entropy
    shannon_entropy = -1*np.sum([p_i*np.log2(p_i) for p_i in freqs if p_i > 0])

    return expected_length, shannon_entropy

def produce_huffman_report(symbols_and_freqs, dictionary):
    """Saves report of Huffman Coding as `report_Huffman.txt`.
    
        If the file exists, the function overwrites. Be careful while calling if the former
        report is required!
    """
    # merge codes with symbols and frequency lise
    symbols_and_freqs_codes = [elem + [dictionary[elem[0]]] for elem in symbols_and_freqs]
    expected_length, shannon_entropy = calculate_expected_lengh_and_shannon_entropy(symbols_and_freqs, dictionary)
    index = 1
    with open(os.path.join('.', 'sample_output', 'report_Huffman.txt'), 'w') as file:
        file.write('Huffman Coding Result Report:\n')
        file.write('-----------------------------\n')
        file.write('Index\tSymbol\tFrequency\tLength\tCode\n')
        for (symbol, freq, code) in symbols_and_freqs_codes:
            code_length = len(code)
            file.write(f'{index}\t{symbol}\t{freq}\t{code_length}\t{code}\n')
            index += 1
        file.write('-----------------------------\n')
        file.write(f'Expected Length = {expected_length}\n')
        file.write(f'Entropy = {shannon_entropy}\n')

def huffman_coded_string_report(text, encoded_text, symbols_and_freqs, dictionary):
    """Generates a file that contains Huffman Coding results and input/encoded strings"""
    symbols_and_freqs_codes = [elem + [dictionary[elem[0]]] for elem in symbols_and_freqs]
    number_of_symbols = len(symbols_and_freqs_codes)
    length_of_encoded_text = len(encoded_text)
    # total_code_length = sum([len(elem) for elem in dictionary.values()])
    expected_length, shannon_entropy = calculate_expected_lengh_and_shannon_entropy(symbols_and_freqs, dictionary)
    with open(os.path.join('.', 'sample_output', 'report_Huffman_coded_string.txt'), 'w') as file:
        file.write('Huffman Coding Result Report:\n')
        file.write('-----------------------------\n')
        file.write(f'Total number of symbols: {number_of_symbols}\n')
        file.write(f'Total number of output bits: {length_of_encoded_text} \n')
        file.write(f'Average symbol length: {expected_length}\n')
        file.write('-----------------------------\n')
        file.write(f'Input String:\n{text}\n')
        file.write('-----------------------------\n')
        file.write(f'Output String:\n{encoded_text}\n')

def huffman_decoded_string_report(encoded_text, decoded_text, symbols_and_freqs, dictionary):
    """Generates a file that contains Huffman Decoding results and encoded/decoded strings"""
    symbols_and_freqs_codes = [elem + [dictionary[elem[0]]] for elem in symbols_and_freqs]
    number_of_symbols = len(symbols_and_freqs_codes)
    length_of_encoded_text = len(encoded_text)
    # total_code_length = sum([len(elem) for elem in dictionary.values()])
    expected_length, shannon_entropy = calculate_expected_lengh_and_shannon_entropy(symbols_and_freqs, dictionary)
    with open(os.path.join('.', 'sample_output', 'report_Huffman_decoded_string.txt'), 'w') as file:
        file.write('Huffman Decoding Result Report:\n')
        file.write('-----------------------------\n')
        file.write(f'Total number of symbols: {number_of_symbols}\n')
        file.write(f'Total number of decoded bits: {length_of_encoded_text} \n')
        file.write(f'Average symbol length: {expected_length}\n')
        file.write('-----------------------------\n')
        file.write(f'Encoded String:\n{encoded_text}\n')
        file.write('-----------------------------\n')
        file.write(f'Decoded String:\n{decoded_text}\n')
