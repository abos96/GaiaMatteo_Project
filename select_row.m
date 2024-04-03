function [r,label] = select_row(txt)
    rws = txt.nRows:-1:1;   % Row number array

    input = inputdlg("Enter row number (max: "+txt.nRows+"):",'Input row',[1 50]);
    input = str2num(input{:});

    if isempty(input) || input<=0 || input>txt.nRows    % Input value must be numeric, positive and less than the number of rows
        errordlg('Incorrect input!','Input Error');
        error('Incorrect input!');
    else
        r = rws(input);
        label = "Row "+int2str(input);
    end
end