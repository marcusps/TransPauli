#! /usr/bin/env python2.5

import csv
import string

def string2pauli( str ):
    """
    This function converts an ASCII string into a multi-qubit Pauli
    operator.
    """
    M = { '1' : 0, 'X' : 1, 'Z' : 2, 'Y' : 3 }
    P = { 'i' : 1, '-' : 2, '-i' : 3 }
    if str[0:1] == '-i':
        p = 3
        str = str[2:]
    elif str[0] == '-':
        p = 2
        str = str[1:]
    elif str[0] == 'i':
        p = 1
        str = str[1:]
    elif str[0] == '+':
        p = 0
        str = str[1:]
    else :
        p = 0
    ilist = map(lambda x: M[x], str)
    return PauliOp( ilist + [ p ] )

class PauliOp:
    
    def __init__( self, tensor_prod_and_sign ):
        # add error checking
        self.plist = tensor_prod_and_sign
        self.N = len( self.plist ) - 1

    def string( self ):
        """ 
        Translate Pauli operator to a string
        """
        M = [ '1', 'X', 'Z', 'Y' ]
        phase = ['+','i','-','-i']
        str = phase[ self.plist[-1] ] + ''.join(map( lambda p: M[p], self.plist[0:-1] ))
        return str

    def cnot( self, ctrl, targ ):
        """
        Apply a controlled-NOT operation
        """
        _p_table = {  0 : [0,0],  1 : [0,1],
                      2 : [2,2],  3 : [2,3],
                      4 : [1,1],  5 : [1,0],
                      6 : [3,3],  7 : [3,2],
                      8 : [2,0],  9 : [2,1],
                     10 : [0,2], 11 : [0,3],
                     12 : [3,1], 13 : [3,0],
                     14 : [1,3], 15 : [1,2] }

        _s_table = {  0 : 0,  1 : 0,
                      2 : 0,  3 : 0,
                      4 : 0,  5 : 0,
                      6 : 2,  7 : 2,
                      8 : 0,  9 : 0,
                     10 : 0, 11 : 0,
                     12 : 0, 13 : 0,
                     14 : 2, 15 : 2 }

        index = (self.plist[ctrl] << 2) ^ self.plist[targ]
        self.plist[-1] = ( self.plist[-1] + _s_table[ index ] ) % 4 
        self.plist[ctrl], self.plist[targ] = _p_table[ index ]

    def csign( self, ctrl, targ ):
        """
        Apply a controlled sign operation
        """
        self.h( targ )
        self.cnot( ctrl, targ )
        self.h( targ )

    def cy( self, ctrl, targ ):
        """
        Apply a controlled Y operation
        """
        self.p( targ )
        self.cnot( ctrl, targ )
        self.p( targ )
        self.p( targ )
        self.p( targ )

    def swap( self, a, b ):
        """
        Swap the operators between qubits a and b
        """
        if a >= self.N or b >= self.N:
            raise 'Index out of bounds'
        self.plist[a], self.plist[b] = self.plist[b], self.plist[a]

    def x( self, a ):
        """
        Apply the Pauli x operator to qubit a
        """
        # flip sign if qubit is Y or Z
        if self.plist[a] & 2:
            self.plist[-1] = (self.plist[-1] + 2) % 4 

    def z( self, a ):
        """
        Apply the Pauli z operator to qubit a
        """
        # flip sign if qubit is X or Y
        if self.plist[a] & 1:
            self.plist[-1] = (self.plist[-1] + 2) % 4 

    def y( self, a ):
        """
        Apply the Pauli y operator to qubit a
        """        
        # flip sign if qubit is X or Z
        if bool( self.plist[a] & 1 ) ^ bool( self.plist[a] & 2 ) :
            self.plist[-1] = (self.plist[-1] + 2) % 4 
        
    def h( self, a ):
        """
        Apply the Hadamard operation to qubit a
        """
        # swap between X and Z, leave Y unchanged
        self.plist[a] = ( self.plist[a] & 1 ) << 1 ^ ( self.plist[a] & 2 ) >> 1
        # flip sign if qubit a is Y
        if self.plist[a] == 3 :
            self.plist[-1] = (self.plist[-1] + 2) % 4

    def p( self, a):
        """
        Apply the Phase operation, diag(1,i), to qubit a 
        """
        self.plist[a] = ( self.plist[a] & 2 ) ^ ( ( self.plist[a] & 1 ) << 1 ) ^ ( self.plist[a] & 1 )
        if self.plist[a] == 1:
            self.plist[-1] = (self.plist[-1] + 2) % 4
            

    def execute_qasm( self, filename ):
        """
        Apply the QASM instructions in the file to the Pauli operator
        (indexing in QASM is 1 based)
        """
        trace = False
        def tron():
            trace = True
        def troff():
            trace = False
        funcs = { 'x' : self.x,
                  'y' : self.y,
                  'z' : self.z,
                  'h' : self.h,
                  'p' : self.p,
                  'cnot' : self.cnot,
                  'cx' : self.cnot,
                  'csign' : self.csign,
                  'cz' : self.csign,
                  'cy' : self.cy,
                  'swap' : self.swap,
                  'tron' : tron,
                  'troff' : troff
                  }
        file = csv.reader(open(filename), delimiter=',')
        row_counter = 0
        for row in file:
            row_counter += 1
            if len(row) > 0 and len(row[0]) > 0:
                # handle inline comments            
                row[-1] = row[-1].split('#')[0]
                # convert all arguments to integers
                arg = map( lambda x: int(x) - 1, row[1:] )
            elif row == [''] or len(row) == 0:
                # do nothing
                continue
            else:
                # handle inline comments
                row[-1] = row[-1].split('#')[0]
                arg = []
            row = map( lambda x: x.strip(), row )
            # now we handle the different commands
            if row[0] == 'tron':
                trace = True
            elif row[0] == 'troff':
                trace = False
            elif len(row[0]) == 0:
                # do nothing
                continue
            else:
                cmd, arg0 = row[0].split()
                arg = [ int(arg0) - 1] + arg 
                if trace:
                    print self.string(), ' ', cmd, '\t', \
                        arg, '\t', \
                        ( funcs[ cmd ]( *arg ), self.string() )[1]
                else:
                    funcs[ cmd ]( *arg )

if __name__ == '__main__':
    """
    Main execution function.

    It will read a Pauli operator from the command line,
    as well as the path of a QASM file describing the 
    circuit/unitary  which acts on the Pauli operator.

    The initial operator can be thought of as a stabilizer
    operator of the input, and we are computing how this 
    stabilizer operator evolves.
    """
    p = string2pauli( '+1Y1XZ' )
    print p.string()
    p.execute_qasm('five.qasm')

# TODO: 
# - add index boundary checks
# - do truth table tests to ensure correctness

