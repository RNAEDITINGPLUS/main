<?php
namespace Mining\Controller;
use Think\Controller;
class QueryController extends Controller {
	public function index(){
		$this->display();
	}
	
	public function verify(){
    	$verify = new \Think\Verify();
    	$verify->useImgBg = true;
    	$verify->entry(1);
    }
	
	public function whats($mail, $job){
		$job = M('jobs');
    	$ts = time();
    	if ($job->autoCheckToken($_POST) && check_verify(I('verify'))) {
			$map['user'] = $mail;
			$map['vocde'] = $job;
			$ret = $job->where($map)->find();
			if ($ret) {
				if ($ret['status'] == 0){
					$this->success('Your job is now waiting for confirmation, please check your mail box and follow the instruction.', '/Mining/Query/index', 5);
				}else if ($ret['status'] == 5){
					Header("Location: https://www.rnaeditplus.org/Viewer/Job/showJob/job/{$ret['vcode']}"); 
				}else{
					$this->success('Your job is running now, this may take some time. And we will mail you, when the task is finished.', '/Mining/Query/index', 5);
				}
			}else{
				$this->error('No matched job.');
			}
		}else{
			$this->error("> Are you attacking RNA Editing Plus server?<br />This form is undes Bro Li's protection.");
		}
	}
}
?>