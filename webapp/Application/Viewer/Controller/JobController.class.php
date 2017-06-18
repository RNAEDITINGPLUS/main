<?php
namespace Viewer\Controller;
use Think\Controller;
class JobController extends Controller {
    public function showJob($job){
    	session("search", '');
		$jobDB = M('jobs');
		$map = array(
			'status'	=>	'5',
			'vcode'		=>	$job,
		);
		
		$job = $jobDB->where($map)->find();
		if ($job){
			$miDB = M('mirediting');
			$u3DB = M('utr3editing');
			$spDB = M('splicing_event');
			$mioDB = M('mir_overview');
			$misDB = M('missense');
			
			/*
			$miDB = M('smir');
			$u3DB = M('sutr3');
			$spDB = M('sse');
			$mioDB = M('mir_overview');
			$misDB = M('smis');
			*/
			$aeDB = M('adar_exp');
			$this->assign('hit', $job);
			session('jobId', $job['trace']);
			$cd = $this->countCandidate($job['trace']);
			$mi = $miDB->where('job='.$job['trace'])->getField('mirna', true);
			$u3 = $u3DB->where('job='.$job['trace'])->getField('gene', true);
			$sp = $spDB->where('job='.$job['trace'])->select();
			$mis = $misDB->where('job='.$job['trace'])->select();
			$ae = $aeDB->where('tissue='.$job['tissue'])->select();
			$mio = $mioDB->where('job_id='.$job['trace'])->select();
			if ($mio){
				$mio = $mio[0];
				session('reg', $mio['reg']);
				session('rel', $mio['rel']);
			}
			$mi = array_unique($mi);
			$u3 = array_unique($u3);
			$this->assign('adar_sites', $ae);
			$this->assign('u3', $u3);
			$this->assign('mi', $mi);
			$this->assign('sp', $sp);
			$this->assign('cd', $cd);
			$this->assign('mio', $mio);
			$this->assign('mis', $mis);
			$this->display();
		}else{
			$this->error('Cannot find your job, please check your input!', U('Mining/Query'));
		}
	}
	
	private function countCandidate($trace){
		$ret = array(
			'sites'	=>	0,
			'genes'	=>	0,
		);
		$sites = M('candidates');
		$map = array(
			'job'	=>	$trace,
		);
		$cds = $sites->where("job='".$trace."'")->getField('refseq_gid', true);
		$ret['sites'] = count($cds);
		$ret['genes'] = count(array_unique($cds));
		return $ret;
	}
	
	private function geneCluster($trace){
		$sites = M('candidates');
		$map = array(
			'job'	=>	$trace,
		);
		$ret = $sites->where($map)->group('RefSeq_gid')->select();
	}
	
	public function saveRecords($job){
		$jobs = M('jobs');
		if ($jobs->autoCheckToken($_POST)){
			$map = array(
				'vcode'	=>	$job,
			);
			$jobId = $jobs->where($map)->getField('trace');
			$Model = new \Think\Model();
			$fileName = create_guid();
			set_time_limit(0);
			$query = "SELECT `region`, `position`, `strand`, `refseq_gid`, `refseq_feat`, `repmask_gid`, `frequency`  FROM  `candidates` WHERE `job`='{$jobId}' INTO OUTFILE 'D://webapp//Outer//{$fileName}.txt';";
			$ret = $Model->execute($query);
			Header("Location: http://www.rnaeditplus.org/Outer/{$fileName}.txt"); 
		}else{
			$this->error("> Are you attacking RNA Editing Plus server?<br />This form is under Bro Li's protection.");
		}
	}
}
?>