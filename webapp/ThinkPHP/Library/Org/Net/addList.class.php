<?php
namespace Org\Net;

class addList{
	public $args0;
	public $args1;
	public $args2;
	public $args3;
	function __construct($inputIds, $idType, $listName, $listType) {
		$this->args0 = $inputIds;
		$this->args1 = $idType;
		$this->args2 = $listName;
		$this->args3 = intval($listType);
	}
}
?>