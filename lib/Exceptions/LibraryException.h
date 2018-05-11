#include <exception> 

class LibraryException : public std::exception 
{
    public: 
        ~LibraryException() throw() {}; 
        const char* what() const throw() {return "An exception regarding the library 'hybridsim' occurred"; } 
};

class ForceException : public LibraryException 
{
        public: 
        ~ForceException() throw() {}; 
        const char* what() const noexcept {return "An overflow error occurred in a force calculation. \n"; }     

}; 

class FENEException : public ForceException 
{
    public: 
        int IndexFirst; 
        int IndexSecond; 
        double Force; 
        std::string ErrorMessage; 
        ~FENEException() throw() {};  
        const char* what() const noexcept {return ErrorMessage.c_str(); }
        FENEException(int First, int Second, double aForce): 
            IndexFirst{First}, 
            IndexSecond {Second}, 
            Force(aForce) {
                ErrorMessage = "An overflow error occurred in a FENE force calculation between the following particles: \n Particle 1: "+std::to_string(IndexFirst)+" \n Particle 2: "+std::to_string(IndexSecond)+" \n Force: "+std::to_string(Force)+"\n" ; 
            }
};    

class RLJException : public ForceException 
{
    public: 
        std::string ErrorMessage; 
        ~RLJException() throw() {};  
        const char* what() const noexcept {return ErrorMessage.c_str(); }
        RLJException(int FirstId, Vector3d FirstPos, Vector3d FirstVel, int SecondId, Vector3d SecondPos, Vector3d SecondVel, double aForce) 
        {
            ErrorMessage = "An overflow error occurred in a Lennard-Jones force calculation between the following particles: \n Particle 1: "+std::to_string(FirstId)+" Position: "+std::to_string(FirstPos(0))+" "+std::to_string(FirstPos(1))+" "+std::to_string(FirstPos(2))+" Velocity: "+std::to_string(FirstVel(0))+" "+std::to_string(FirstVel(1))+ " "+std::to_string(FirstVel(2))+" \n Particle 2: "+std::to_string(SecondId)+" "+std::to_string(SecondPos(0))+" "+std::to_string(SecondPos(1))+" "+std::to_string(SecondPos(2))+" Velocity: "+std::to_string(SecondVel(0))+" "+std::to_string(SecondVel(1))+ " "+std::to_string(SecondVel(2))+" \n Force: "+std::to_string(aForce)+"\n" ; 
        }
}; 

class LJSurfException : public ForceException 
{
    public: 
        std::string ErrorMessage; 
        ~LJSurfException() throw() {};  
        const char* what() const noexcept {return ErrorMessage.c_str(); }
        LJSurfException(int Id, Vector3d Pos, Vector3d Vel, double aForce) 
        {
         ErrorMessage = "An overflow error occurred in a Lennard-Jones surface force calculation for the following particle: "+std::to_string(Id)+" Position: "+std::to_string(Pos(0))+" "+std::to_string(Pos(1))+" "+std::to_string(Pos(2))+" Velocity: "+std::to_string(Vel(0))+" "+std::to_string(Vel(1))+ " "+std::to_string(Vel(2))+" \n Force: "+std::to_string(aForce)+"\n" ; 
        }
};

class CellAllocationException : public std::exception {
    public: 
        std::string ErrorMessage; 
        ~CellAllocationException() throw() {}; 
        const char* what() const noexcept {return ErrorMessage.c_str(); }
        CellAllocationException(const Particle part, const std::array<int, 3> CellIndex) {
            ErrorMessage = "Bad allocation of Particle "+std::to_string(part.Identifier)+" with Position {"+std::to_string(part.Position(0))+", "+std::to_string(part.Position(1))+", "+std::to_string(part.Position(2))+"} at CellIndex {"+std::to_string(CellIndex[0])+", "+std::to_string(CellIndex[1])+", "+std::to_string(CellIndex[2])+"}"; 
        } 
}; 

