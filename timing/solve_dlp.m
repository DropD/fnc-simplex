function [x, fval, exitflag, output] = solve_dlp(fname)
    lp = fileread(fname);
    fst = regexp(lp, 'Maximize') + length('Maximize') + 1;
    fen = regexp(lp, 'Subject') - 2;
    abst = fen + length('Subject To') + 3;
    aben = regexp(lp, 'End') - 2;

    f = -str2num(lp(fst:fen));
    ab = str2num(regexprep(lp(abst:aben), '<= ', ''));

    A = ab(:,1:end-1);
    b = ab(:,end);

    [x, fval, exitflag, ouput] = linprog(f, A, b, [], [], [], [], [], optimset('Display', 'iter'));
end
