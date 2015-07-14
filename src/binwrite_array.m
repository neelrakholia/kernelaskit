function fid=binwrite_array( data, filename)
% data : double array
% filename 
% you should call fclose(fid)

endiancheck_string = '1234ABCD';
endiancheck = hex2dec(endiancheck_string);

fid = fopen( filename, 'w');
fwrite(fid, endiancheck, 'int64');
fwrite(fid, data(:), 'double');




