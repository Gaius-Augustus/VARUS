#ifndef VNAME_H
#define VNAME_H

#ifndef VNAME

#define VNAME(x) #x;

#endif // VNAME

#ifndef VDUMP

#define VDUMP(x) std::cout << #x << ":\t" << x << std::endl

#endif // VDUMP

#endif//VNAME
