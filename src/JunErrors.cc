#include "JunErrors.hh"

void MiaoError(string errorStr)
{
  cout<<YELLOW<<BOLD<<"# Miao! "<<COLOR_END;
  cout<<RED<<BOLD<<"Error: "<<COLOR_END;
  cout<<errorStr<<endl;
  exit(0);
}

void MiaoWarn(string warnStr)
{
  cout<<YELLOW<<BOLD<<"# Miao! "<<COLOR_END;
  cout<<BLUE<<BOLD<<"Warn: "<<COLOR_END;
  cout<<warnStr<<endl;
}

void MiaoHint(string hintStr)
{
  cout<<YELLOW<<BOLD<<"# Miao! "<<COLOR_END;
  cout<<GREEN<<BOLD<<"Hint: "<<COLOR_END;
  cout<<hintStr<<endl;
}