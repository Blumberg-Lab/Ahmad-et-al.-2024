function txtNUM = jkAddZero2Num(digitSpan, inputNum)
%To add 0s to txtNum for easy numbering
%Code by Jangjin Kim, 2012-08-28
%
%Input param definition
%digitSpan [1 x 1]              how many zeros will be added to txt number
%inputNum [1 x 1]               input number to be converted to text
%
%Output param definition
%txtNUM [1 x 1]                 text converted number with 0s as fillers

numZEROS = digitSpan - size(num2str(inputNum), 2);
txtZEROS = repmat(['0'], 1, numZEROS);

txtNUM = [txtZEROS num2str(inputNum)];