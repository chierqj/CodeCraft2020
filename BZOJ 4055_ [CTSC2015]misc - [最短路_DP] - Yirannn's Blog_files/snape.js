(function() {

	'use strict';

	$(".Share").click(function() {
		$(".share-wrap").fadeToggle("slow");
	});

	$('.qrcode').each(function(index, el) {
		var url = $(this).data('url');
		if ($.fn.qrcode) {
			$(this).qrcode({
				text: url,
				width: 150,
				height: 150,
			});
		}
	});

	$(document).on("click", ".Love",
		function() {
			$(this).postLike();
	});
	$.fn.postLike = function() {
		if ($(this).hasClass('done')) {
			alert('您已赞过该文章');
			return false;
		} else {
			$(this).addClass('done');
			var id = $(this).data("id"),
			action = $(this).data('action'),
			rateHolder = $(this).children('.count');
			var ajax_data = {
				action: "love",
				um_id: id,
				um_action: action
			};
			$.post("/wp-admin/admin-ajax.php", ajax_data,
			function(data) {
				$(rateHolder).html(data);
			});
			return false;
		}
	};

	var offset = 300,
		offset_opacity = 1200,
		scroll_top_duration = 700,
		$back_to_top = $('.cd-top');
	$(window).scroll(function() {
		($(this).scrollTop() > offset) ? $back_to_top.addClass('cd-is-visible') : $back_to_top.removeClass('cd-is-visible cd-fade-out');
		if ($(this).scrollTop() > offset_opacity) {
			$back_to_top.addClass('cd-fade-out');
		}
	});

	$back_to_top.on('click', function(event) {
		event.preventDefault();
		$('body,html').animate({
			scrollTop: 0,
		}, scroll_top_duration);
	});

	console.log("╔╦╗┬ ┬┌─┐┌┬┐┌─┐  ╔═╗┌┐┌┌─┐┌─┐┌─┐  ╔╦╗┌─┐┌┬┐┌─┐  ╔╗ ┬ ┬  ╦  ╦┌┬┐┬─┐┌─┐┬┌─┐\n ║ ├─┤├┤ │││├┤   ╚═╗│││├─┤├─┘├┤   ║║║├─┤ ││├┤   ╠╩╗└┬┘  ╚╗╔╝ │ ├┬┘│ ││└─┐\n ╩ ┴ ┴└─┘┴ ┴└─┘  ╚═╝┘└┘┴ ┴┴  └─┘  ╩ ╩┴ ┴─┴┘└─┘  ╚═╝ ┴    ╚╝  ┴ ┴└─└─┘┴└─┘");
	console.log("Snape 主题下载：https://github.com/Vtrois/Snape");
	console.log("Snape 主题使用：https://www.vtrois.com/projects/snape-faq.html");
	console.log("Snape 文章样式：https://www.vtrois.com/projects/snape-article-style.html");
}());