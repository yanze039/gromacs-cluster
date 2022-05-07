import os 
import glob 
import json
import shutil
import numpy as np
import pandas as pd
import argparse
# import MDAnalysis as mda
# from MDAnalysis.analysis import rms

from pathlib import Path
# from dpdispatcher import Machine, Resources, Task, Submission
import subprocess

def get_name(name):
    name_idx = 0
    while os.path.exists("{}.{:03d}".format(name, name_idx)):
        name_idx += 1
    return "{}.{:03d}".format(name, name_idx)

def movedir(src, end):
    if not os.path.exists(end):
        os.mkdir(end)
    file_list = glob.glob(os.path.join(src, "*"))
    for ff in file_list:
        shutil.move(ff, end)

def clean_dpdispatcher(clean_path):
    all_task = []
    tag_list = ["*.json", "*_task_fail", "*.out", "*_finished", "*.sub", "*_job_id", "dpdispatcher.log"]
    for tag in tag_list:
        all_task += glob.glob(os.path.join(clean_path, tag))
    for ff in all_task:
        os.remove(ff)

class Executer:
    def __init__(self):
        pass
        
    def execute(self, cmd, stdin=None, work_path="."):
        subp = subprocess.Popen(
            args=cmd,
            stdin =subprocess.PIPE,
            cwd=work_path.absolute(),
            # stdout=subprocess.PIPE,
            # stderr=subprocess.PIPE,
            encoding="utf-8",
            shell=True
        )
        if stdin is not None:
            subp.stdin.write(stdin)
            subp.stdin.close()
        subp.wait()
        assert subp.returncode == 0
        return subp.returncode

class Gromacs(Executer):
    def __init__(self, config):
        super(Gromacs, self).__init__()
        self.__dict__.update(config)
    
    def get_group_type(self, cluster_group):
        if cluster_group == 'protein':
            return 1
        elif cluster_group == 'protein-h':
            return 2
        elif cluster_group == 'calpha':
            return 3
        elif cluster_group == 'backbone':
            return 4
        else:
            raise RuntimeError("Unsupported Group Type.")

    def cluster(self, pdb_list):
        tmp_dir_idx = 0
        while(os.path.exists("tmp.{:03d}".format(tmp_dir_idx))):
            tmp_dir_idx += 1
        tmp_dir = "tmp.{:03d}".format(tmp_dir_idx)
        self.tmp_dir = tmp_dir
        os.mkdir(tmp_dir)
        allpdb = open(f"{tmp_dir}/all.pdb", "w")
        for pdb in pdb_list:
            with open(pdb, "r") as content:
                allpdb.write(content.read())
            allpdb.write("\n")
        allpdb.close()
        pdb2trr_cmd = "gmx trjconv -f all.pdb -s {} -o traj.trr -timestep 1".format(pdb_list[0])
        cls_cmd = "gmx cluster -f traj.trr -s {} -method {} -cutoff {} -sz -clid -cl".format(
                    pdb_list[0], self.cluster_method, self.cluster_cutoff
                )
        cluster_group_idx = self.get_group_type(self.cluster_group)
        self.execute(pdb2trr_cmd, stdin="1\n", work_path=Path(tmp_dir))
        self.execute(cls_cmd, stdin="{}\n1\n".format(cluster_group_idx), work_path=Path(tmp_dir))
        return cls_cmd
    
    # def em(self, pdb_list, work_dir, task_config):
    #     work_dir = os.path.join(work_dir, "relaxed_model")
    #     if os.path.exists(work_dir):
    #         output_bak = get_name(work_dir)
    #         movedir(work_dir, output_bak)
    #         print("Previous results are backed up at {}".format(output_bak))
    #         shutil.rmtree(work_dir)
    #     restrain_group = self.get_group_type(self.restrain_group)
    #     os.mkdir(work_dir)
    #     task_list = []
    #     tmp_dir_list = []
    #     for pdb in pdb_list:
    #         p = Path(pdb)
    #         sub_dir = os.path.join(work_dir, p.stem)
    #         os.mkdir(sub_dir)
    #         shutil.copyfile("mdp/minim.mdp", os.path.join(sub_dir, "minim.mdp"))
    #         # shutil.copytree(self.forcefield, os.path.join(sub_dir, os.path.basename(self.forcefield)))
    #         cmd =  'echo "5\n1\n" | gmx pdb2gmx -f {} -o processed.gro '.format(pdb)
    #         cmd += '&& echo "{}\n" | gmx genrestr -f processed.gro -o posre.itp -fc {} {} {} '.format(restrain_group, self.kappa, self.kappa, self.kappa)
    #         cmd += "&& gmx grompp -f minim.mdp -c processed.gro -p topol.top -r processed.gro -o em.tpr -maxwarn 2 "
    #         cmd += "&& gmx mdrun -deffnm em -v -ntmpi 1 -nt 4 "
    #         cmd += '&& echo "1\n" | gmx trjconv -s em.tpr -f em.gro -o relaxed_{}.pdb -pbc mol '.format(p.stem)
    #         cmd += '&& mv relaxed_{}.pdb ../'.format(p.stem)
    #         print(cmd)
    #         task_list.append(
    #             Task(command=cmd, task_work_path=p.stem)
    #         )
    #         tmp_dir_list.append(os.path.join(work_dir, p.stem))
    #     resources = Resources.load_from_dict(task_config["resources"])
    #     machine = Machine.load_from_dict(task_config["machine"])
    #     submission = Submission(
    #         work_base=work_dir,
    #         machine=machine, 
    #         resources=resources,
    #         task_list=task_list, 
    #     )
    #     submission.run_submission()
    #     clean_dpdispatcher(work_dir)
    #     if self.clean:
    #         for tmp in tmp_dir_list:
    #             shutil.rmtree(tmp)
    #     return
    
    def read_cls_log(self):
        cls_info = {}
        with open(os.path.join(self.tmp_dir, "cluster.log"), "r") as clslog:
            current_cls = None
            for line in clslog.readlines():
                if (not "|" in line) or ("cl." in line):
                    continue
                else:
                    info = [xx.split() for xx in line.split("|")]
                    if len(info[0]) > 0:
                        cls_info[info[0][0]] = {
                            "population": int(info[1][0]),
                            "cluster_center": int(info[2][0]),
                            "cluster_member": [int(xx) for xx in info[3]]
                        }
                        current_cls = info[0][0]
                    else:
                        cls_info[current_cls]["cluster_member"] += [int(xx) for xx in info[3]]
        if self.clean:
            shutil.rmtree(self.tmp_dir)
        return cls_info

def gather_json(json_list):
    all_data = {}
    for jj in json_list:
        with open(jj, "r") as jd:
            jdata = json.load(jd)
        all_data.update(jdata)
    return all_data

# def cal_rmsd(ref, model, group="backbone"):
#     u = mda.Universe(model) 
#     ref = mda.Universe(ref) 
#     rmsd_value = rms.rmsd(u.select_atoms(group).positions,  # coordinates to align
#                     ref.select_atoms(group).positions,  # reference coordinates
#                     center=True,  # subtract the center of geometry
#                     superposition=True)  # superimpose coordinates
#     return rmsd_value * 0.1

def main():
    parser = argparse.ArgumentParser(description="Clustering pdb files within a path")
    parser.add_argument("config", type=str, default="config.json",
                        help="Config file for clustering in JSON format.")
    parser.add_argument("--output", "-o", type=str, default="cluster_result", dest="output",
                        help="The output directory.")
    parser.add_argument("--prefix", "-p", type=str, default="model", dest="prefix",
                        help="The prefix of collected pdb files. (model/multimer) ")
    parser.add_argument("--score", "-s", type=str, default="plddt", dest="score",
                        help="The score type used to rank classes. (ptm/plddt)")
    parser.add_argument("--test", action='store_true', default=False, dest="test",
                        help="Wether to use test mode.")
    # parser.add_argument("--unclean", action='store_true', default=False, dest="unclean",
    #                     help="wether to reserve the tmp files.")
    args = parser.parse_args()
    
    with open(args.config, "r") as cfg:
        jdata = json.load(cfg)

    ## Clustering 
    gmx = Gromacs(jdata["gromacs"])
    pdb_list = glob.glob(
        os.path.join(jdata["inference_dir"], "{}*.pdb".format(args.prefix))
    )
    if args.test:
        pdb_list = pdb_list[:50]
    else:
        pdb_list = pdb_list
    pdb_list.sort()
    if len(pdb_list) == 0:
        raise RuntimeError(f"No pdb file matches {args.prefix}*.pdb format at path {args.inference_dir}!")
    mapping = {}
    for idx, pdb in enumerate(pdb_list):
        mapping[idx] = pdb
    cls_cmd = gmx.cluster(pdb_list)
    cls_info = gmx.read_cls_log()
    
    ## Gathering score info
    json_list = glob.glob(
        os.path.join(jdata["inference_dir"], f"*{args.score}.json")
    )
    if len(json_list) == 0:
        print(f"No *{args.score}.json file exists within dir {args.inference_dir}")
        raise RuntimeError(f"Please check scoreing JSON at path {args.inference_dir}.")
    all_score = gather_json(json_list)

    # make result dir
    if os.path.exists(args.output):
        output_bak = get_name(args.output)
        movedir(args.output, output_bak)
        print("Previous results are backed up at {}".format(output_bak))
        shutil.rmtree(args.output)
    os.mkdir(args.output)
    
    ## Process All Information
    ## Generate results
    final_info = {}
    cls_msg = "{} clusters are found under cutoff {} nm".format(len(cls_info.keys()), gmx.cluster_cutoff)
    csv_cls_info = {}
    
    paired_cls_lddt_list = []
    for cluster in cls_info.keys():
        _info = {}
        _info2 = {}

        mapped_center = os.path.basename(mapping[cls_info[cluster]["cluster_center"]])
        mapped_member = [
            os.path.basename(mapping[xx]) for xx in cls_info[cluster]["cluster_member"]
        ]

        mapped_member_with_score = [
            (xx, float(all_score[xx])) for xx in mapped_member
        ]
        ranked_score_list = sorted(mapped_member_with_score, key=lambda item:item[1], reverse=True)
        score_list = [float(all_score[xx]) for xx in mapped_member]
        
        _info["population"] = cls_info[cluster]["population"]
        _info[f"max {args.score} in cls"] = ranked_score_list[0]
        _info[f"min {args.score} in cls"] = ranked_score_list[-1]
        _info[f"mean {args.score} in cls"] = np.mean(score_list)
        _info[f"mid {args.score} in cls"] = np.median(score_list)
        _info[f"cluster_center"] = (mapped_center, float(all_score[mapped_center]))
        _info[f"cluster_member"] = ranked_score_list
        _info2[f"population"] = cls_info[cluster]["population"]
        _info2[f"max {args.score} in cls"] = np.max(score_list)
        _info2[f"min {args.score} in cls"] = np.min(score_list)
        _info2[f"mean {args.score} in cls"] = np.mean(score_list)
        _info2[f"mid {args.score} in cls"] = np.median(score_list)
        # rmsd_range = 0
        paired_cls_lddt_list.append((cluster, np.max(score_list)))
        
        final_info[cluster] = _info
        csv_cls_info[cluster] = _info2
    
    ranked_paired_cls_lddt_list = sorted(paired_cls_lddt_list, key=lambda item:item[1], reverse=True)
    ranked_final_info = {}
    ranked_csv_cls_info = {}
    ranked_final_info["summary"] = cls_msg
    for ii in range(len(ranked_paired_cls_lddt_list)):
        tag = f"Rank.{args.score}.{ii}"
        ranked_final_info[tag] = final_info[ranked_paired_cls_lddt_list[ii][0]]
        ranked_csv_cls_info[tag] = csv_cls_info[ranked_paired_cls_lddt_list[ii][0]]

    ranked_final_info["cluster_cmd"] = cls_cmd
    print("Clustering Done.")
    
    copy_tag = ["best", 'worst', 'center']
    for cluster in ranked_csv_cls_info:
        subcluster_dir = os.path.join(args.output, cluster)
        os.mkdir(subcluster_dir)
        _info = ranked_final_info[cluster]
        for ii in range(int(np.min([len(copy_tag), int(_info["population"])]))):
            if copy_tag[ii] == "best":
                shutil.copyfile(
                        os.path.join(jdata["inference_dir"], _info[f"max {args.score} in cls"][0]),
                        os.path.join(subcluster_dir, "best_"+_info[f"max {args.score} in cls"][0])
                    )
            elif copy_tag[ii] == "worst":
                shutil.copyfile(
                        os.path.join(jdata["inference_dir"], _info[f"min {args.score} in cls"][0]),
                        os.path.join(subcluster_dir, "worst_"+_info[f"min {args.score} in cls"][0])
                    )
            elif copy_tag[ii] == "center":
                shutil.copyfile(
                        os.path.join(jdata["inference_dir"], _info["cluster_center"][0]),
                        os.path.join(subcluster_dir, "center_"+_info["cluster_center"][0])
                    )


    # output JSON
    with open(os.path.join(args.output, "cluster_result.json"), "w") as clsres:
        json.dump(ranked_final_info, clsres, indent=4)
    print(cls_msg)
    
    # output csv
    df = pd.DataFrame(ranked_csv_cls_info)
    df.T.to_csv(os.path.join(args.output, "cluster_result.csv"))

    
if __name__ == "__main__":
    main()
    