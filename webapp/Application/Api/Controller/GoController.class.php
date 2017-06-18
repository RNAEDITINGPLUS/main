<?php
namespace Api\Controller;
use Think\Controller;
class GoController extends Controller {
	public function mirRedirction($mir){
		$queDB = M('mirna_map');
		$map = array('mir' => $mir, );
		$ret = $queDB->where($map)->find();
		if ($ret){
			redirect('http://www.mirbase.org/cgi-bin/mature.pl?mature_acc='.$ret['acc']);
		}else{
			$this->error('No miRNA record!');
		}
	}
}
?>