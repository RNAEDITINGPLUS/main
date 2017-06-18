<?php
function enhencedGeneSearch($term){
	$geneInfo = M('gene_info');
	$map['synonyms'] = array('LIKE', array($term.'%', '%|'.$term.'|%', '%|'.$term), 'OR');
	$data = $geneInfo->where($map)->field('symbol, synonyms')->select();
	return $data;
}

function translator4Splicing($t){
	$dict = array(
		3	=>	"inactived (or weakened) 3' splicing site",
		5	=>	"inactived (or weakened) 5' splicing site",
		6	=>	"enhanced 3' splicing site",
		9	=>	"weakened splicing",
		10	=>	"enhanced 5' splicing site",
		12	=>	"new 3' splicing site",
		15	=>	"weakened 5' splicing site",
		16	=>	"new branch site",
		17	=>	"inactived branch site",
		18	=>	"weakened branch site",
		19	=>	"enhenced branch site",
		20	=>	"new 5' splicing",
	);
	return $dict[$t];
}

function translator4Seed($is_seed){
	$icon = "";
	if ($is_seed){
		$icon = '<i class="glyphicon glyphicon-ok"></i>';
	}else{
	$icon = '<i class="glyphicon glyphicon-remove"></i>';
	}
	return $icon;
}