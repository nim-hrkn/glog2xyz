import sys
import argparse
import os

class scf_result:
    def __init__(self, input_orientation_atoms,standard_orientation_atoms, 
         original_axis_forces, energy, maxcycle ):
         self.input_orientation_atoms = input_orientation_atoms 
         self.standard_orientation_atoms = standard_orientation_atoms
         self.original_axis_forces = original_axis_forces
         self.energy = energy
         self.maxcycle = maxcycle
    def xyz1(self):
         """
         convert to xyz format
         use input_orientation_atoms and original_axis_forces
         """
         s = []
         n = len(self.input_orientation_atoms)
         s.append(str(n))
         s.append("input_orientation_atoms energy %s maxcycle %s" %(self.energy,self.maxcycle))
         for x,y in zip(self.input_orientation_atoms,self.original_axis_forces):
             flag = False
             if x[0]!=y[0]:
                 flag = True
             if x[1]!=y[1]:
                 flag = True
             if flag:
                 print "internal error, atom index inconsistent"
                 print x
                 print y
                 raise Exception
             x.extend(y[2:])
             s.append( " ".join(x))
         return "\n".join(s)

    def xyz2(self):
         """
         convert to xyz format
         use standard_orientation_atoms, no force for them
         """
         s = []
         n = len(self.standard_orientation_atoms)
         s.append(str(n))
         s.append("standard_orientation_atoms energy %s maxcycle %s" %(self.energy,self.maxcycle))
         for x in self.standard_orientation_atoms:
             s.append( " ".join(x))
         return "\n".join(s)

    def xyz3(self):
         """
         convert to xyz format
         use input_orientation_atoms, no force for them
         """
         s = []
         n = len(self.input_orientation_atoms)
         s.append(str(n))
         s.append("input_orientation_atoms energy %s maxcycle %s" %(self.energy,self.maxcycle))
         for x in self.input_orientation_atoms:
             s.append( " ".join(x))
         return "\n".join(s)


    def xyz(self):
         self.xyztype = 0 
         #print self.original_axis_forces is None, self.standard_orientation_atoms is None, self.input_orientation_atoms is None

         if self.original_axis_forces is not None:
             if self.input_orientation_atoms is not None:
                 self.xyztype = 1
                 return self.xyz1()
             if self.standard_orientation_atoms is not None:
                 self.xyztype = 2
                 return self.xyz2()
         else:
             if self.standard_orientation_atoms is not None:
                 self.xyztype = 2
                 return self.xyz2()
             elif self.input_orientation_atoms is not None:
                 self.xyztype = 3
                 return self.xyz3()
             else:
                 print "xyz type error"
                 raise Exception

    def __str__(self):
         """
         string 
         """
         s = []
         s.append("input orientation")
         for x,y in zip(self.input_orientation_atoms,self.original_axis_forces):
             flag = False
             if x[0]!=y[0]:
                 flag = True
             if x[1]!=y[1]:
                 flag = True
             if flag:
                 print "internal error, atom index inconsistent"
                 print x
                 print y
                 raise Exception
             x.extend(y[2:])
             s.append( " ".join(x))
         if False:
           s.append("standard orientation")
           for x in self.standard_orientation_atoms:
             s.append( " ".join(x))
         #s.append("original axis forces")
         #for x in self.original_axis_forces:
         #    s.append( " ".join(x))
         s.append("energy = %s" % self.energy)
         s.append("maxcycle = %s" % self.maxcycle)
         return "\n".join(s)

class glog: 
    def __init__(self):
        """
        @param filename: log file name
        """
        arg = self.argparse()
        self.filename = arg.logfile
        self.xyzfile = arg.xyz 

        self.natom = None
        lines = self.load_log(self.filename)
        self.status = self.check_normal_exit(lines)
        if self.status:
            self.process_all(lines)

    def load_log(self,filename):
        """
        real gauusian log file
        @param filename: log filename of gaussian
        @return : log text
        """
        with open(filename,"rb") as f:
            data = f.read()
        lines = data.split("\n")
        return lines

    def process_atoms(self,it):
        """
        process atomic coodinates
        @param it: iterator of lines
        @return : it , atoms , it:terator, atoms: list of atomic positions
        """
        key1 = " -------"
        n_key1 = len(key1)
        if True:
            x = it.next()
            x = it.next()
            x = it.next()
            x = it.next()
        atoms = []
        if self.natom is None:
            while True:
                x = it.next()
                if x[:n_key1]== key1:
                    break 
                else:
                    atoms.append(x.split())
        else:
           for i in range(self.natom):
               x = it.next()
               atoms.append(x.split())
        return it, atoms 
        
    def process_scf_done(self,line):
        """
        process SCF Done
        @param line: lines of the files
        @return : done,energy, maxcycles , done=True|False, energy=total energy, maxcycles=cycles of iterations
        """
        done = False
        maxcycles = 0 
        s = line.split()
        i = s.index("A.U.")
        energy = s[i-1]
        if s[i+1] == "after" and s[i+3] == "cycles":
            maxcycles = s[i+2] 
            done = True
        return done,energy,maxcycles

    def check_normal_exit(self,lines):
        """
        check Normal termination
        @param lines: lines of the file
        @return : flag = True|False
        """
        key4 = " Normal termination of"
        n_key4 = len(key4)
        #print "number of lines = ",len(lines)
        x = lines[-2]
        #print "last line (",x,")"
        if x[:n_key4 ] == key4:
                return True
        return False

    def process_forces(self,it):
        """
        process force lines
        @param it: iterator
        @return :  a list of forces
        """
        x = it.next()
        x = it.next()
        key1 = " ---"
        n_key1 = len(key1)
        forces = []
        while True:
            x = it.next()
            if x[:n_key1] == key1:
                break
            forces.append(x.split())
        return forces

    def check_force_lines(self,lines):
        key5 = " ***** Axes restored to original set"
        n_key5 = len(key5)
        for x in lines:
             if x[:n_key5] == key5:
                 self.force_lines = True
                 return True
        self.force_lines = False
        return False

    def process_one_scf_force(self,it):
        """
        process log file
        @param lines: log file
        @return: none
        """
        key1 = "                          Input orientation:"
        n_key1 = len(key1)
        key2 = "                         Standard orientation:"
        n_key2 = len(key2)
        key3 = " SCF Done:"
        n_key3 = len(key3)
        key4 = " Normal termination of"
        n_key4 = len(key4)
        key5 = " ***** Axes restored to original set"
        n_key5 = len(key5)
        key6 = " Center     Atomic                   Forces"
        n_key6 = len(key6)
      
        input_orientation_atoms = None
        standard_orientation_atoms = None
        original_axis_forces = None
        while True:
            try:
               x = it.next()
            except:
               break 
            if x[:n_key1] == key1:
                it, input_orientation_atoms = self.process_atoms(it) 
            elif x[:n_key2 ] == key2: 
                it, standard_orientation_atoms = self.process_atoms(it) 
            elif x[:n_key3 ] == key3: 
                flag, energy, maxcycle = self.process_scf_done(x)
            elif x[:n_key4 ] == key4: 
                return False, it, None
            elif x[:n_key5 ] == key5: 
                x = it.next()
                x = it.next()
                if x[:n_key6] == key6:
                    original_axis_forces = self.process_forces(it)
                    return flag,it, scf_result(input_orientation_atoms, standard_orientation_atoms,original_axis_forces,energy, maxcycle)

        return False,it,None

    def process_one_scf_noforce(self,it):
        """
        process log file
        @param lines: log file
        @return: none
        """
        key1 = "                          Input orientation:"
        n_key1 = len(key1)
        key2 = "                         Standard orientation:"
        n_key2 = len(key2)
        key3 = " SCF Done:"
        n_key3 = len(key3)
        key4 = " Normal termination of"
        n_key4 = len(key4)
      
        input_orientation_atoms = None
        standard_orientation_atoms = None
        original_axis_forces = None
        while True:
            try:
               x = it.next()
            except:
               break 
            if x[:n_key1] == key1:
                it, input_orientation_atoms = self.process_atoms(it) 
            elif x[:n_key2 ] == key2: 
                it, standard_orientation_atoms = self.process_atoms(it) 
            elif x[:n_key3 ] == key3: 
                flag, energy, maxcycle = self.process_scf_done(x)
                return flag,it, scf_result(input_orientation_atoms, standard_orientation_atoms,original_axis_forces,energy, maxcycle)

            elif x[:n_key4 ] == key4: 
                return False, it, None

        return False,it,None


    def process_all(self,lines):
        """
        process all the prodecures
        @param lines: lines of the file
        """
        self.atoms = []
        if self.check_force_lines(lines):
            it = iter(lines)
            while True:
                flag, it, status = self.process_one_scf_force(it) 
                if not flag: 
                    break
                self.atoms.append( status )
        else:
            it = iter(lines)
            while True:
                flag, it, status = self.process_one_scf_noforce(it)
                if not flag:
                    break
                self.atoms.append( status )


    def xyz(self):
        """
        output xyz format to the file
        @return : xyz formatted text
        """
        if not self.status:
            return 
        lines = []
        for a in self.atoms :
             y = a.xyz()
             xyztype = a.xyztype
             lines.append(y)
        s = "\n".join(lines)+"\n"
        with open(self.xyzfile,"wb") as f:
            f.write(s)
        self.xyztype = xyztype
        print "output a xyz file to ", self.xyzfile

    def make_xyzfilename(self,logpath):
        p,ext = os.path.splitext(logpath)        
        return p+".xyz"

    def argparse(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-xyz",help="xyz file name (output)",default=None)  
        parser.add_argument("logfile", help="log file of Gaussian")  
        args = parser.parse_args()

        if args.xyz == None:
            args.xyz = self.make_xyzfilename(args.logfile)

        return args

if __name__ == "__main__":
 
    g = glog()
    g.xyz()

    if g.status:
         returncode = 0
         if g.xyztype !=1:
             print "without forces"
             returncode = 1
         sys.exit( returncode )
    else:
         print "Not Normal Termnation"
         sys.exit( -1 )


