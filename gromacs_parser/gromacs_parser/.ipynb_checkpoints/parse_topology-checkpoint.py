import os
import subprocess

class ParseTopology:
    
    def __init__(self,topology_path):
        self.topology_path=topology_path
        self.gmx_top_dir =self.get_gmx_top_dir()
        self.topology_text=self.read_topology(self.topology_path)
        
    def get_gmx_top_dir(self):
        stdout=subprocess.run("which gmx",capture_output=True,text=True,shell=True).stdout
        if stdout!="":
            gmx_path=stdout.replace("\n","")
            gmx_lib_dir="/".join(stdout.split("/")[:-2])+'/share/gromacs/top'
            if os.path.isdir(gmx_lib_dir):
                return gmx_lib_dir
            else:
                print("Warning: gromacs top dir not found please assign self.gmx_top_dir")
                return ""
        else:
            print("Warning: gromacs top dir not found please assign self.gmx_top_dir")
            return ""
        
    def trim_data(self, data):
        #trimming comments
        no_comments_data = [line.split(";")[0] for line in data.split("\n")]
        #get rid of empty lines
        no_empty_lines = [ line for line in no_comments_data if line.split()!=[]]
        trimmed_data="\n".join(no_empty_lines)
        return trimmed_data

    def read_recuresive(self, dir_path, file_name):
        data_lines=[]
        file_path=dir_path+"/"+file_name
        with open(file_path,'r') as h:
            data=self.trim_data(h.read())
        for line in data.split("\n"):
            if "#include" in line:
                include_file=line.split('"')[1]
                include_path=dir_path+'/'+include_file
                if not os.path.isfile(include_path):
                    include_path=self.gmx_top_dir+'/'+include_file
                include_dir="/".join(include_path.split('/')[:-1])
                include_file_name=include_path.split("/")[-1]
                include_lines=self.read_recuresive(include_dir,include_file_name)
                for include_line in include_lines:
                    data_lines.append(include_line)
            else:
                data_lines.append(line)
        return data_lines

    def read_topology(self, topology_path):
        topology_dir="/".join(topology_path.split('/')[:-1])
        topology_file=topology_path.split("/")[-1]
        whole_data="\n".join(self.read_recuresive(topology_dir, topology_file))
        return whole_data



    
