function vals = Select_Rows(f, i)
    
    S    = size(f,1);
    vals = f((i-1)*S + (1:S)');
    
end