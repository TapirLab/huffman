"""
Function that reads the text file and extract the relative frequencies

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

function [sorted_frequencies, symbol_list_descending_order] = ReadTextInput(FileName)
    fid = fopen(FileName); % Initialize the file with the pointer
    B = textscan(fid, '%s'); % Read the file and dump into the variable
    fclose(fid); % Close the file
    textinput = strjoin(B{1}); % Join the strings into a horizontal vector
    unique_list_of_input = unique(textinput); % Extract the unique list of symbols
    for(k = 1 : length(unique_list_of_input)) % Iterate through the unique symbol list
        relative_frequency(k) = numel(find(unique_list_of_input(k) == textinput))/length(textinput); % Compute the relative frequencies of each symbol
    end
    
    [sorted_frequencies, sorting_idx] = sort(relative_frequency, 'descend'); % Extract the indices for sorting
    symbol_list_descending_order = unique_list_of_input(sorting_idx); % Sort out the unique symbols in descending order

end

