function [] = LeapYear(year)
% LeapYear Determines if the given year is a leap year
    if (mod(year,400)==0 || mod(year,4)==0 && mod(year,100)~=0)
        disp([num2str(year), ' is a leap year']);
    else
        disp([num2str(year), ' is not a leap year']);
    end
end