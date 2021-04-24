"""
Matlab implementation of Huffman Coding.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Huffman Coding
%% -------------------
%% $Author: Tapir Lab.$,
%% $Date: January 5th, 2020$,
%% $Revision: 1.0$
%% $Tapir Lab.$
%% $Copyright: Tapir Lab.$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""

clear all; % Clear the workspace
close all; % Close all open objects
clc; % Clear the terminal

%%%%%%% Definitions %%%%%%%
TreeType = 1; % Defines the weight of the left leaf of any node.
symbol_start_ASCII = 65; % Default symbols are selected to be ASCII characters
maximum_number_of_symbols = 26; % Maximum number of symbols allowed
DecayRate = 0.01; % Decay rate for the exponential distribution
NumberOfRealizations = 50; % In case symbol frequencies are generated in this module, define the number of realizations
number_of_symbols = 6; % Number of symbols to be generated in case random text is produced by this module
TypeOfDistribution = 0; % Default type of distribution in case symbol frequencies are generated in this module; [0]
Scramble = 1; % Decide whether generated symbols are scrambled (set for scamble; reset for no scramble)
alphabet_letters = char([symbol_start_ASCII : symbol_start_ASCII + maximum_number_of_symbols - 1]'); % Form the symbols to be Huffman-coded in capital Latin letters
%%%% Definitions - END %%%%

relative_frequency_list_in_descending_order = GenerateRandomSymbolFrequencies(NumberOfRealizations,TypeOfDistribution,number_of_symbols,DecayRate); % Generate random symbols acconrding to the parameters provided in the preamble
out = GenerateRandomSymbolsAccordingToPredefinedFrequencies(relative_frequency_list_in_descending_order, alphabet_letters, NumberOfRealizations, Scramble);
%%%%%%%%%%%%%% Precomputed Variables %%%%%%%%%%%%%%
number_of_letters = length(relative_frequency_list_in_descending_order); % Identify the number of relative frequencies provided/generated
alphabet = [1 : number_of_letters]'; % Create index for the symbols to keep track of moves and states
%%%%%%%%%%%% Precomputed Variables - END %%%%%%%%%%

%%%%%%%%%%%%%% Initialization List %%%%%%%%%%%%%%
list_of_moves_of_leaves = relative_frequency_list_in_descending_order(:); % Initialize the list of relative frequencies
leaf_node_index_list = [alphabet]; % Keep the list of leaves with appropriate indexing
relative_leaf_position_list = zeros(number_of_letters,number_of_letters); % Initialize the relative leaf positions being left or right
k = 0; % Initialize the state counter
[left, right] = DecideTreeType(TreeType); % Assign the weights for left and right node
%%%%%%%%%%% Initialization List - END %%%%%%%%%%%




input_message = GenerateRandomSymbolsAccordingToPredefinedFrequencies(relative_frequency_list_in_descending_order, alphabet_letters, NumberOfRealizations, Scramble);

% relative_frequency_list_in_descending_order = [0.20
% 0.18
% 0.10
% 0.10
% 0.10
% 0.06
% 0.06
% 0.04
% 0.04
% 0.04
% 0.04
% 0.03
% 0.01];







while(length(list_of_moves_of_leaves) > 1) % Exhaust the frequency list until probability/relative frequencies sum up to unity
    k = k + 1; % Increment the counter that stores the current state
    [lowest_element, second_lowest_element, sum_of_lowest_frequencies] = PickTheTwoLowestElements(list_of_moves_of_leaves); % Pick the two lowest items in the frequency list and obtain their symbol indices (original indices) along with their sum for further states
    [larger_freq_idx, insertion_point_idx, lower_freq_idx] = IdentifyTheSplitLists(list_of_moves_of_leaves, sum_of_lowest_frequencies); % Identify how many sublists form when sum of the two lowest items is obtained
        if(isempty(larger_freq_idx)) % In case the sum of the two lowest frequencies exceeds the maximum frequency value in the contemporary state (no split in the list occurs)
            if(sum_of_lowest_frequencies == 1) % Check whether the final state is reached. In case it is, execute the following; if not skip on ELSE
                relative_leaf_position_list = PopulateTreeTable(leaf_node_index_list, 0, relative_leaf_position_list, k, 1, left); % Populate the indices of the contributors that will be represented on the left leaf/node. The second parameter should be reset (binary 0) since all indices present (except for the contributors to the sum) will be shifted one unit down as a block 
                relative_leaf_position_list = PopulateTreeTable(leaf_node_index_list, 0, relative_leaf_position_list, k, 2, right); % Populate the indices of the contributors that will be represented on the right leaf/node. The second parameter should be reset (binary 0) since all indices present (except for the contributors to the sum) will be shifted one unit down as a block 
                [list_of_moves_of_leaves] = CleanUp(sum_of_lowest_frequencies,list_of_moves_of_leaves); % Reduce the list with appropriate indexing
            else % The sum of the two lowest frequencies exceeds the maximum frequency value in the contemporary list but final state is not reached (no split in the list occurs, all of the indices except for the contributors to the two lowest will be shifted one unit down)
                list_of_indices_to_be_shifted = []; % Allocate the memory for the indices to be shifted one unit down
                [relative_leaf_position_list, leaf_node_index_list, list_of_indices_to_be_shifted_from_below] = PopulateTreeTable(leaf_node_index_list, 1, relative_leaf_position_list, k, lowest_element, left); % Populate the indices of the contributors that will be represented on the left leaf/node. The second parameter should be set (binary 1) since all indices present (except for the contributors to the sum) will be shifted one unit down as a block 
                [relative_leaf_position_list, leaf_node_index_list, list_of_indices_to_be_shifted_from_above] = PopulateTreeTable(leaf_node_index_list, 1, relative_leaf_position_list, k, second_lowest_element, right); % Populate the indices of the contributors that will be represented on the right leaf/node.  The second parameter should be set (binary 1) since all indices present (except for the contributors to the sum) will be shifted one unit down as a block
                [leaf_node_index_list] = UpdateListOfTreeTable(leaf_node_index_list,alphabet, list_of_indices_to_be_shifted_from_below, list_of_indices_to_be_shifted_from_above, k); % Update the states based on the information above
                [list_of_moves_of_leaves] = CleanUp([sum_of_lowest_frequencies; list_of_moves_of_leaves(1 : second_lowest_element - 1)],list_of_moves_of_leaves); % Reduce the list with appropriate indexing
            end
        else
            larger_freq_idx_original = FindIndexWithModularFixForZero(leaf_node_index_list, k, larger_freq_idx, number_of_letters); % Identify the indices for the items that remain above the split point      
            lower_freq_idx_original = FindIndexWithModularFixForZero(leaf_node_index_list, k, lower_freq_idx, number_of_letters); % Identify the indices for the items that remain below the split point        
            [relative_leaf_position_list, tmp1, tmp2, new_idx_final_original] = PopulateTreeTable(leaf_node_index_list, 0, relative_leaf_position_list, k, lowest_element, left); % Populate the indices of the contributors that will be represented on the left leaf/node. The second parameter should be reset (binary 0) since all indices present (except for the contributors to the sum) will be shifted one unit down as a block  
            [relative_leaf_position_list, tmp1, tmp2, new_idx_final_previous_original] = PopulateTreeTable(leaf_node_index_list, 0, relative_leaf_position_list, k, second_lowest_element, right); % Populate the indices of the contributors that will be represented on the right leaf/node. The second parameter should be reset (binary 0) since all indices present (except for the contributors to the sum) will be shifted one unit down as a block  
            merge = SplitAndMergeList(list_of_moves_of_leaves(lower_freq_idx), list_of_moves_of_leaves(larger_freq_idx), sum_of_lowest_frequencies); % Split the current list of frequency list and form the new list by inserting the sum of two lowest frequencies in the current state
            [leaf_node_index_list] = AssignNewStates(leaf_node_index_list, k, larger_freq_idx_original, insertion_point_idx, lower_freq_idx_original, new_idx_final_previous_original,new_idx_final_original); % Reorganize the indices to keep track of the evolution of all leaves/nodes
            [list_of_moves_of_leaves] = CleanUp(merge,list_of_moves_of_leaves); % Reduce the list with appropriate indexing
        end    
    end
    
    

code_dictionary = AssignBinaryCode(relative_leaf_position_list, right, alphabet_letters); % Assigning binary codes to leaves based on the table entries obtained throughout the process
[entropy, report_list] = EvaluateResults(code_dictionary, relative_frequency_list_in_descending_order, alphabet_letters); % Generate the report for the dictionary and input alphabet 

coded_message = ApplyHuffmanCompression(input_message, code_dictionary, alphabet_letters(alphabet)) % Apply Huffman to encode the input string with the dictionary given

decoded_string = HuffmanDecoder(coded_message, code_dictionary, alphabet_letters(alphabet));

disp(['Coding complete. Check the report file.']);


%%%% Function block that computes the entropy, forms report containing results, and creates a text file to store them all %%%%  
function [entropy, report] = EvaluateResults(Dictionary, FrequencyList, AlphabetSymbols)
report = []; % Allocate slot for the report table for workspace (This line and corresponding return paramter could be omitted)
file_pointer = fopen('report_Huffman.txt','w'); % Initialize the file with the pointer
fprintf(file_pointer,'Huffman Coding Result Report:\n'); % Write the title of the report
fprintf(file_pointer,'--------\n'); % Write the separator in the file for visual enhancement
fprintf(file_pointer,'%s\t %s\t %s\t %s\t\t %s\n','Index', 'Symbol', 'Frequency','Length', 'Code'); % Write the headings for the pattern

    for(k = 1 : length(FrequencyList)) % Iterate through the symbols provided
        report(k,1) = FrequencyList(k); % Obtain the k-th frequency value to write into the file
        report(k,2) = Dictionary(k).strlength; % Obtain the k-th binary code to be assigned to k-th symbol
        fprintf(file_pointer,'%d\t %s\t %1.2f\t\t %d\t\t %s\n', k, AlphabetSymbols(k), report(k,1), report(k,2), Dictionary(k)); % Write entries into the file according to the pattern provided
    end
    entropy = report(:,1)'*report(:,2); % Calculate the entropy
    fprintf(file_pointer,'--------\n'); % Write the separator in the file for visual enhancement
    fprintf(file_pointer,'Entropy: %1.4f bit \n',entropy); % Write the computed entropy value in the file
fclose(file_pointer); % Close the file accessed by the pointer

end
            
%%%% Function block that populates the indices in appropriate places for the leaf-node evolution throughout Huffman coding structure %%%%  
function [TreeTable, Ordered_List_Positions, Index_Shift, Indices] = PopulateTreeTable(Ordered_List_Positions, Ordered_List_Position_Update, TreeTable, Current_State_Index, Search_Value, Leaf_Position_Value)
    Indices = find(Search_Value == Ordered_List_Positions(:, Current_State_Index)); % Find the index of the value to be sought for for the current state
    Next_State_Index = Current_State_Index + 1; % Identify the next state via index
    Index_Shift = []; % Allocate the slot for the indices to be shifted
    for(uu = 1 : length(Indices)) % Iterate through the number of indices found
        if(TreeTable(Indices(uu), Next_State_Index) == 0) % In case zero entries found, execute the following
            TreeTable(Indices(uu), Next_State_Index) = Leaf_Position_Value; % Place the value of the leaf (Left or Right for binary tree) in the correct index
            if(Ordered_List_Position_Update) % In case the sum of the lowest two frequency items reach the highest frequency (this is detected by setting the second parameter of this function call)
                Ordered_List_Positions(Indices(uu), Next_State_Index) = Ordered_List_Position_Update; % Place the correspnding entry on top by setting (assigning "1") the next state index
            end
            Index_Shift = [Index_Shift ; Indices(uu)]; % Accumulate the indices to be processed for further use
        end
    end  
end

%%%% Function block that cleans up the frequency list by reducing its size in light of the recently merged list %%%%
function [output_array] = CleanUp(input_array, output_array)
    tmp_array = [input_array]; % Assign temporary variable
    output_array = []; % Reallocate an empty slot for the output array
    output_array = tmp_array; % Re-assign the input with clean input (in MATLAB size reduction is required, otherwise output_array might lead to inaccessible indices)
end


%%%% Function block that finds the indices %%%%
function [output_list] = FindIndexWithModularFixForZero(Ordered_List_Positions, Current_State_Index, Search_Value, size_of_alphabet)
    Search_Value = Search_Value(:); % Vectorize the input that contains the indices to be sought for. Make sure column vector is obtained. It will be row vector in the upcoming lines
    Search_Space = Ordered_List_Positions(:, Current_State_Index); % Construct the search space in column vector format
    output_list = mod(find(Search_Value' == Search_Space(:)),size_of_alphabet); % Search row vector in the column vector. The output is given by a column vector whose entries are in modular form. So, in MATLAB, find([1 3 4] == [1 2 3 4 5]') yields [1 8 14]'.
    zeroth_idx = find(output_list == 0); % Identify the entries (which will be used as indices) with zero values. Due to modular structure, entries delivered with zero should be identified. This line could be replaced by a FOR-LOOP.
    if(~isempty(zeroth_idx)) % In case zero valued entry (which will be used as indices) is detected, execute the following line
        output_list(zeroth_idx) = size_of_alphabet; % Replace the zero valued entries with the total number of symbols fed in
    end
end


%%%% Function block that reorganizes the indices regarding the frequency update %%%%
function [output] = SplitAndMergeList(Lower_List, Upper_List, Middle_List)
    split_up = Upper_List(:); % Vectorize the items that store the list of indices that represent the symbols contributing to higher frequency values
    split_down = Lower_List(:); % Vectorize the items that store the list of indices that represent the symbols contributing to lower frequency values
    middle_list = Middle_List(:); % Vectorize the items that store the list of indices that represent the symbols equal to the sum of the two lowest frequency values
    output = [split_up; middle_list; split_down]; % Construct the list of indices based on the frequency update
end


%%%% Function block that picks the two lowest frequency items in the frequency list %%%%
function [last_item, last_item_above, sum_of_freqs] = PickTheTwoLowestElements(FrequencyList)
    last_item = length(FrequencyList); % The lowest frequency index is the length of the list of the contemporary frequency states
    last_item_above = last_item - 1; % The second lowest frequency index
    sum_of_freqs = FrequencyList(last_item) + FrequencyList(last_item_above); % Sum of the two lowest frequencies
end

%%%% Function block that split the frequency list based on the previous state in the light of the sum of the two lowest frequencies %%%%
function [ListOfAboveItems, SplitItemIndex, ListOfBelowItems] = IdentifyTheSplitLists(FrequencyList, Threshold)
    ListOfAboveItems = find(Threshold <= FrequencyList); % Identify the indices that are greater than or equal to the threshold (sum of the two lowest frequencies in the previous state)
    SplitItemIndex = max(ListOfAboveItems) + 1; % Locate the index where the sum of the two lowest frequency will be placed in the next state
    ListOfBelowItems = [SplitItemIndex : length(FrequencyList) - 2]'; %Identify the indices that are lower than the sum of the two lowest frequencies in the previous state by excluding the indices of the two lowest frequencies in the previous state
end

%%%% Function block that shifts the symbol indices one unit down when sum of the two lowest frequencies becomes the highest frequency in the list %%%%
function [ListOfUpdatedIndices] = UpdateListOfTreeTable(ListOfUpdatedIndices,ListOfAlphabetSymbols, ListOfIndicesToBeUpdatedBelow, ListOfIndicesToBeUpdatedAbove, CurrentState)
    NextState = CurrentState + 1; % Define the index for the next state
    ListOfIndexToBeUpdated = setdiff(ListOfAlphabetSymbols,[ListOfIndicesToBeUpdatedBelow(:) ; ListOfIndicesToBeUpdatedAbove(:)]); % Identify the indices that will be shifted down in the list of updated frequencies
    ListOfUpdatedIndices(ListOfIndexToBeUpdated, NextState) = ListOfUpdatedIndices(ListOfIndexToBeUpdated, CurrentState) + 1; % Shift the indices one unit down
end

%%%% Function block that assigns binary code to each symbol %%%%            
function dictionary = AssignBinaryCode(TreeTable, Leaf_Node_Binary_Value, AlphabetSymbols)
dictionary = ""; % Initialize the string that will store the bit sequences for each symbol decided by Huffman coding
  for(k = 1 : size(TreeTable,2)) % Iterate through the table that holds for the tree constructed by Huffman coding row-wise where k-th row holds all evolution of k-th symbol through each and every update process
    [vals, idx]= find(TreeTable(k, 1 : end) > 0); % Identify the indices that show activity for the k-th symbol (no activity entry is represented with zero)
    bit_sequence_in_reverse_order = (TreeTable(k, idx) == Leaf_Node_Binary_Value); % Mask the activity-demonstrating entries with [Leaf_Node_Binary_Value]: matching entries yield binary 1, non-matching entries yield binary 0. Note that iterations (k) start from 1 implying tree is scanned through in reverse order. 
    bit_sequence = num2str(fliplr(bit_sequence_in_reverse_order)); % Binary code is flipped around to reflect the prefix code places on most significant bits and converted into strings for easy representation
    bit_sequence(bit_sequence == ' ') = []; % Remove the white spaces that will appear after applying [num2str()] function
    dictionary(k) = string(bit_sequence); % Assign the binary code for the k-th symbol
  end
  
end

%%%% Function block that assign new states (positions) for each symbol upon updating the list with combining two lowest frequencies %%%%            
function [LeafTrackList] = AssignNewStates(LeafTrackList, CurrentState, ListOfAboveItems, SplitItemIndex, ListOfBelowItems, IndexOfOriginalSecondLowest, IndexOfOriginalLowest)
    NextState = CurrentState + 1; % Define the index for the next state 
    LeafTrackList(ListOfBelowItems, NextState) = LeafTrackList(ListOfBelowItems,CurrentState) + 1; % Shift the items that have lower frequencies one unit down in the updated list of frequencies due to insertion of super node
    LeafTrackList(ListOfAboveItems, NextState) = LeafTrackList(ListOfAboveItems, CurrentState); % Keep the items that have higher frequencies with indices unchanged
    LeafTrackList(IndexOfOriginalLowest, NextState) = SplitItemIndex; % Assign the new index to position the symbol on the splitting point of the list of frequencies (Symbol that has relatively higher frequency)
    LeafTrackList(IndexOfOriginalSecondLowest, NextState) = SplitItemIndex; % Assign the new index to position the symbol on the splitting point of the list of frequencies (Symbol that has relatively lower frequency)
end


%%%% Function block that generates the relative frequencies of the symbols %%%%
function [RelativeFrequencies] = GenerateRandomSymbolFrequencies(NumberOfRealizations, DistributionType, NumberOfSymbols, DecayRate)
    switch DistributionType % Select the distribution type of the random variable
        case 0 % Uniform distribution
            rv = randi(NumberOfSymbols, NumberOfRealizations,1); % Realize a set of uniformly distributed random samples    
        case 1 % Exponential distribution
            rv = exprnd(DecayRate, NumberOfRealizations, 1); % Realize a set of exponentially distributed random samples
        % case 2 % Fill out the block in case some other distributions are
        % included in here
    end    
    [freqs,bins] = hist(rv(:), NumberOfSymbols); % Form the histogram with the number of bins equal to the number of symbols in the dictionary
    RelativeFrequencies = sort(freqs./sum(freqs),'descend'); % Provide relative frequencies in descending order for Huffman coder
end

%%%% Function block that decides the weights of the leaves %%%%
function [WeightOfLeftNode, WeightOfRightNode] = DecideTreeType(TreeType) 
    if(TreeType > 0) % Decide whether the leaves are right-handed or left-handed
        WeightOfRightNode = 10; % Right node has higher weight
        WeightOfLeftNode = 1; % Set a nonnegative constant for left node
    else
        WeightOfRightNode = 1; % % Right node has lower weight
        WeightOfLeftNode = 10; % Set a nonnegative constant for left node
    end
end


function out = GenerateRandomSymbolsAccordingToPredefinedFrequencies(FrequencyList, Symbols, NumberOfRealizations, Scramble)
output_text = []; % Allocate slot for the output
    for(k = 1 : length(FrequencyList)) % Iterate through the list of frequencies given in descending order
        tmp = repmat(Symbols(k), 1, round(FrequencyList(k)*NumberOfRealizations)); % Concatenate the symbol of interest in light of the corresponding relative frequency
        output_text = [output_text tmp]; % Concatenate the generated symbol with the previouslu generated ones
    end
    if(Scramble) % Decide whether the output list of symbols to be scrambled or not
        out = output_text(randperm(length(output_text))); % Obtain a random permutation of the input string
    end
end


%%%% Function block that converts input symbols (letters) into bit streams via dictionary %%%%
function binary_string_out = ApplyHuffmanCompression(InputMessage, Dictionary, ListOfSymbols)
    binary_string = ''; % Allocate empty string
    for(k = 1 : length(InputMessage)) % Iterate through the input message symbol by symbol
        symbol = InputMessage(k); % Retrieve k-th symbol
        idx_symbol = find(symbol == ListOfSymbols); % Find the index that corresponds to the mapping
        binary_string = [binary_string Dictionary(idx_symbol)]; % Concatenate the binary coded strings back to back
    end
    binary_string_out = strcat(binary_string{1:end}); % Merge the entire string chunks into a single block
    
    file_pointer = fopen('report_Huffman_coded_string.txt','w'); % Initialize the file with the pointer
    fprintf(file_pointer,'Huffman Coding Result Report:\n'); % Write the title of the report
    fprintf(file_pointer,'Total number of symbols:\t%d\n',length(InputMessage)); % Write the total number of symbols
    fprintf(file_pointer,'Total number of output bits:\t%d\n', length(binary_string_out)); % Write the total number of bits
    fprintf(file_pointer,'Average symbol length:\t%2.4f\n', length(binary_string_out)/length(InputMessage)); % Write the average symbol length
    fprintf(file_pointer,'--------\n'); % Write the separator in the file for visual enhancement
    fprintf(file_pointer,'Input String: %s\n', InputMessage); % Write the input string
    fprintf(file_pointer,'--------\n'); % Write the separator in the file for visual enhancement
    fprintf(file_pointer,'Output String: %s\n', binary_string_out); % Write the output string
    fclose(file_pointer); % Close the file accessed by the pointer    
end


%%%% Function block that decodes the input bit stream %%%
function decoded_stream = HuffmanDecoder(InputBinaryStream, Dictionary, Alphabet)
buffer = []; % Allocate the buffer
decoded_stream = []; % Allocate space for decoded stream
maximum_length_of_dictionary_entries = max(Dictionary.strlength); % Identify the maximum depth of the decoding tree
    for(k = 1 : length(InputBinaryStream)) % Iterate through the input bit stream bit by bit
        buffer = [buffer InputBinaryStream(k)]; % Fill the buffer with bit streams from the left
        if(length(buffer) > maximum_length_of_dictionary_entries) % In case buffer length exceeds the depth of the dictionary
            error('Unrecognized prefix. Decoding failed!'); % An unrecognized prefix is encountered. Decoding should stop.
        else % In case buffer length is not exceeded, execute fthe following
            search_string_idx = find(buffer == Dictionary); % Locate index of the string stored in buffer
            if(~isempty(search_string_idx)) % In case an index returns
                decoded_stream = [decoded_stream Alphabet(search_string_idx)]; % Concatenate the decoded bit string into corresponding symbol in dictionary
                buffer = []; % Empty the buffer for the next symbol
            end
        end

    end
    file_pointer = fopen('report_Huffman_decoded_string.txt','w'); % Initialize the file with the pointer
    fprintf(file_pointer,'Huffman Decoding Result Report:\n'); % Write the title of the report
    fprintf(file_pointer,'Total number of symbols:\t%d\n',length(decoded_stream)); % Write the total number of symbols
    fprintf(file_pointer,'Total number of decoded bits:\t%d\n', length(InputBinaryStream)); % Write the total number of bits
    fprintf(file_pointer,'--------\n'); % Write the separator in the file for visual enhancement
    fprintf(file_pointer,'Input String: %s\n', InputBinaryStream); % Write the input string
    fprintf(file_pointer,'--------\n'); % Write the separator in the file for visual enhancement
    fprintf(file_pointer,'Output String: %s\n', decoded_stream); % Write the output string
    fclose(file_pointer); % Close the file accessed by the pointer      

end