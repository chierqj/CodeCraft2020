(function($) {
    "use strict";
    var bootsnav = {
        initialize: function() {
            this.event();
            this.hoverDropdown();
            this.navbarSticky();
            this.navbarScrollspy()
        },
        event: function() {
            var getNav = $("nav.navbar.bootsnav");
            if(getNav.hasClass("brand-center")) {
                var postsArr = new Array(),
                    index = $("nav.brand-center"),
                    $postsList = index.find('ul.navbar-nav');
                index.find('ul.navbar-nav > li').each(function() {
                    postsArr.push($(this).html())
                });
                var firstList = postsArr.splice(0, Math.round(postsArr.length / 2)),
                    secondList = postsArr,
                    ListHTML = '';
                var createHTML = function(list) {
                    ListHTML = '';
                    for(var i = 0; i < list.length; i++) {
                        ListHTML += '<li>' + list[i] + '</li>'
                    }
                }
                createHTML(firstList);
                $postsList.html(ListHTML);
                index.find("ul.nav").first().addClass("navbar-left");
                createHTML(secondList);
                $postsList.after('<ul class="nav navbar-nav"></ul>').next().html(ListHTML);
                index.find("ul.nav").last().addClass("navbar-right");
                index.find("ul.nav.navbar-left").wrap("<div class='col-half left'></div>");
                index.find("ul.nav.navbar-right").wrap("<div class='col-half right'></div>");
                index.find('ul.navbar-nav > li').each(function() {
                    var dropDown = $("ul.dropdown-menu", this),
                        megaMenu = $("ul.megamenu-content", this);
                    dropDown.closest("li").addClass("dropdown");
                    megaMenu.closest("li").addClass("megamenu-fw")
                })
            }
            if(getNav.hasClass("navbar-sidebar")) {
                $("body").addClass("wrap-nav-sidebar");
                getNav.wrapInner("<div class='scroller'></div>")
            } else {
                $(".bootsnav").addClass("on")
            }
            if(getNav.find("ul.nav").hasClass("navbar-center")) {
                getNav.addClass("menu-center")
            }
            if(getNav.hasClass("navbar-full")) {
                $("nav.navbar.bootsnav").find("ul.nav").wrap("<div class='wrap-full-menu'></div>");
                $(".wrap-full-menu").wrap("<div class='nav-full'></div>");
                $("ul.nav.navbar-nav").prepend("<li class='close-full-menu'><a href='#'><i class='fa fa-times'></i></a></li>")
            } else if(getNav.hasClass("navbar-mobile")) {
                getNav.removeClass("no-full")
            } else {
                getNav.addClass("no-full")
            }
            if(getNav.hasClass("navbar-mobile")) {
                $('.navbar-collapse').on('shown.bs.collapse', function() {
                    $("body").addClass("side-right")
                });
                $('.navbar-collapse').on('hide.bs.collapse', function() {
                    $("body").removeClass("side-right")
                });
                $(window).on("resize", function() {
                    $("body").removeClass("side-right")
                })
            }
            if(getNav.hasClass("no-background")) {
                $(window).on("scroll", function() {
                    var scrollTop = $(window).scrollTop();
                    if(scrollTop > 34) {
                        $(".navbar-fixed").removeClass("no-background")
                    } else {
                        $(".navbar-fixed").addClass("no-background")
                    }
                })
            }
            if(getNav.hasClass("navbar-transparent")) {
                $(window).on("scroll", function() {
                    var scrollTop = $(window).scrollTop();
                    if(scrollTop > 34) {
                        $(".navbar-fixed").removeClass("navbar-transparent")
                    } else {
                        $(".navbar-fixed").addClass("navbar-transparent")
                    }
                })
            }
            $(".btn-cart").on("click", function(e) {
                e.stopPropagation()
            });
            $("nav.navbar.bootsnav .attr-nav").each(function() {
                $("li.search > a", this).on("click", function(e) {
                    e.preventDefault();
                    $(".top-search").slideToggle()
                })
            });
            $(".input-group-addon.close-search").on("click", function() {
                $(".top-search").slideUp()
            });
            $("nav.navbar.bootsnav .attr-nav").each(function() {
                $("li.side-menu > a", this).on("click", function(e) {
                    e.preventDefault();
                    $("nav.navbar.bootsnav > .side").toggleClass("on");
                    $("body").toggleClass("on-side")
                })
            });
            $(".side .close-side").on("click", function(e) {
                e.preventDefault();
                $("nav.navbar.bootsnav > .side").removeClass("on");
                $("body").removeClass("on-side")
            });
            $(".tab-nav li").each(function() {
                $("a", this).off("click");
                $("a", this).on("click", function(e) {
                    e.preventDefault();
                    var cek = $(this).closest("li").hasClass("active");
                    if(cek) {
                        $("a", this).off("click")
                    } else {
                        var getID = $(this).attr("href");
                        $(".tab-nav li").removeClass("active");
                        $(this).closest("li").addClass("active")
                    }
                })
            });
            $(".dropdown-tabs .tab-pane").each(function() {
                var getBG = $(this).data("background");
                $(this).css("background-image", "url('" + getBG + "')");
                $(this).addClass("animated")
            });
            $("body").wrapInner("<div class='wrapper'></div>")
        },
        hoverDropdown: function() {
            var getNav = $("nav.navbar.bootsnav"),
                getWindow = $(window).width(),
                getHeight = $(window).height(),
                getIn = getNav.find("ul.nav").data("in"),
                getOut = getNav.find("ul.nav").data("out");
            if(getWindow < 981) {
                $(".scroller").css("height", "auto");
                $("nav.navbar.bootsnav ul.nav").find("li.dropdown").off("mouseenter");
                $("nav.navbar.bootsnav ul.nav").find("li.dropdown").off("mouseleave");
                $("nav.navbar.bootsnav ul.nav").find(".title").off("mouseenter");
                $("nav.navbar.bootsnav ul.nav").off("mouseleave");
                $(".navbar-collapse").removeClass("animated");
                $("nav.navbar.bootsnav ul.nav").each(function() {
                    $(".dropdown-menu", this).addClass("animated");
                    $(".dropdown-menu", this).removeClass(getOut);
                    $("a.dropdown-toggle", this).off('click');
                    $("a.dropdown-toggle", this).on('click', function(e) {
                        e.stopPropagation();
                        $(this).closest("li.dropdown").find(".dropdown-menu").first().stop().fadeToggle().toggleClass(getIn);
                        $(this).closest("li.dropdown").first().toggleClass("on");
                        return false
                    });
                    $('li.dropdown', this).each(function() {
                        $(this).find(".dropdown-menu").stop().fadeOut();
                        $(this).on('hidden.bs.dropdown', function() {
                            $(this).find(".dropdown-menu").stop().fadeOut()
                        });
                        return false
                    });
                    $(".megamenu-fw", this).each(function() {
                        $(".col-menu", this).each(function() {
                            $(".content", this).addClass("animated");
                            $(".content", this).stop().fadeOut();
                            $(".title", this).off("click");
                            $(".title", this).on("click", function() {
                                $(this).closest(".col-menu").find(".content").stop().fadeToggle().addClass(getIn);
                                $(this).closest(".col-menu").toggleClass("on");
                                return false
                            });
                            $(".content", this).on("click", function(e) {
                                e.stopPropagation()
                            })
                        })
                    })
                });
                var cleanOpen = function() {
                    $('li.dropdown', this).removeClass("on");
                    $(".dropdown-menu", this).stop().fadeOut();
                    $(".dropdown-menu", this).removeClass(getIn);
                    $(".col-menu", this).removeClass("on");
                    $(".col-menu .content", this).stop().fadeOut();
                    $(".col-menu .content", this).removeClass(getIn)
                }
                $("nav.navbar.bootsnav").on("mouseleave", function() {
                    cleanOpen()
                });
                $("nav.navbar.bootsnav .attr-nav").each(function() {
                    $(".dropdown-menu", this).removeClass("animated");
                    $("li.dropdown", this).off("mouseenter");
                    $("li.dropdown", this).off("mouseleave");
                    $("a.dropdown-toggle", this).off('click');
                    $("a.dropdown-toggle", this).on('click', function(e) {
                        e.stopPropagation();
                        $(this).closest("li.dropdown").find(".dropdown-menu").first().stop().fadeToggle();
                        $(".navbar-toggle").each(function() {
                            $(".fa", this).removeClass("fa-times");
                            $(".fa", this).addClass("fa-bars");
                            $(".navbar-collapse").removeClass("in");
                            $(".navbar-collapse").removeClass("on")
                        })
                    });
                    $(this).on("mouseleave", function() {
                        $(".dropdown-menu", this).stop().fadeOut();
                        $("li.dropdown", this).removeClass("on");
                        return false
                    })
                });
                $(".navbar-toggle").each(function() {
                    $(this).off("click");
                    $(this).on("click", function() {
                        $(".fa", this).toggleClass("fa-bars");
                        $(".fa", this).toggleClass("fa-times");
                        cleanOpen()
                    })
                })
            } else if(getWindow > 981) {
                $(".scroller").css("height", getHeight + "px");
                if(getNav.hasClass("navbar-sidebar")) {
                    $("nav.navbar.bootsnav ul.nav").each(function() {
                        $("a.dropdown-toggle", this).off('click');
                        $("a.dropdown-toggle", this).on('click', function(e) {
                            e.stopPropagation()
                        });
                        $(".dropdown-menu", this).addClass("animated");
                        $("li.dropdown", this).on("mouseenter", function() {
                            $(".dropdown-menu", this).eq(0).removeClass(getOut);
                            $(".dropdown-menu", this).eq(0).stop().fadeIn().addClass(getIn);
                            $(this).addClass("on");
                            return false
                        });
                        $(".col-menu").each(function() {
                            $(".content", this).addClass("animated");
                            $(".title", this).on("mouseenter", function() {
                                $(this).closest(".col-menu").find(".content").stop().fadeIn().addClass(getIn);
                                $(this).closest(".col-menu").addClass("on");
                                return false
                            })
                        });
                        $(this).on("mouseleave", function() {
                            $(".dropdown-menu", this).stop().removeClass(getIn);
                            $(".dropdown-menu", this).stop().addClass(getOut).fadeOut();
                            $(".col-menu", this).find(".content").stop().fadeOut().removeClass(getIn);
                            $(".col-menu", this).removeClass("on");
                            $("li.dropdown", this).removeClass("on");
                            return false
                        })
                    })
                } else {
                    $("nav.navbar.bootsnav ul.nav").each(function() {
                        $("a.dropdown-toggle", this).off('click');
                        $("a.dropdown-toggle", this).on('click', function(e) {
                            e.stopPropagation()
                        });
                        $(".megamenu-fw", this).each(function() {
                            $(".title", this).off("click");
                            $("a.dropdown-toggle", this).off("click");
                            $(".content").removeClass("animated")
                        });
                        $(".dropdown-menu", this).addClass("animated");
                        $("li.dropdown", this).on("mouseenter", function() {
                            $(".dropdown-menu", this).eq(0).removeClass(getOut);
                            $(".dropdown-menu", this).eq(0).stop().fadeIn().addClass(getIn);
                            $(this).addClass("on");
                            return false
                        });
                        $("li.dropdown", this).on("mouseleave", function() {
                            $(".dropdown-menu", this).eq(0).removeClass(getIn);
                            $(".dropdown-menu", this).eq(0).stop().fadeOut().addClass(getOut);
                            $(this).removeClass("on")
                        });
                        $(this).on("mouseleave", function() {
                            $(".dropdown-menu", this).removeClass(getIn);
                            $(".dropdown-menu", this).eq(0).stop().fadeOut().addClass(getOut);
                            $("li.dropdown", this).removeClass("on");
                            return false
                        })
                    })
                }
                $("nav.navbar.bootsnav .attr-nav").each(function() {
                    $("a.dropdown-toggle", this).off('click');
                    $("a.dropdown-toggle", this).on('click', function(e) {
                        e.stopPropagation()
                    });
                    $(".dropdown-menu", this).addClass("animated");
                    $("li.dropdown", this).on("mouseenter", function() {
                        $(".dropdown-menu", this).eq(0).removeClass(getOut);
                        $(".dropdown-menu", this).eq(0).stop().fadeIn().addClass(getIn);
                        $(this).addClass("on");
                        return false
                    });
                    $("li.dropdown", this).on("mouseleave", function() {
                        $(".dropdown-menu", this).eq(0).removeClass(getIn);
                        $(".dropdown-menu", this).eq(0).stop().fadeOut().addClass(getOut);
                        $(this).removeClass("on")
                    });
                    $(this).on("mouseleave", function() {
                        $(".dropdown-menu", this).removeClass(getIn);
                        $(".dropdown-menu", this).eq(0).stop().fadeOut().addClass(getOut);
                        $("li.dropdown", this).removeClass("on");
                        return false
                    })
                })
            }
            if(getNav.hasClass("navbar-full")) {
                var windowHeight = $(window).height(),
                    windowWidth = $(window).width();
                $(".nav-full").css("height", windowHeight + "px");
                $(".wrap-full-menu").css("height", windowHeight + "px");
                $(".wrap-full-menu").css("width", windowWidth + "px");
                $(".navbar-collapse").addClass("animated");
                $(".navbar-toggle").each(function() {
                    var getId = $(this).data("target");
                    $(this).off("click");
                    $(this).on("click", function(e) {
                        e.preventDefault();
                        $(getId).removeClass(getOut);
                        $(getId).addClass("in");
                        $(getId).addClass(getIn);
                        return false
                    });
                    $("li.close-full-menu").on("click", function(e) {
                        e.preventDefault();
                        $(getId).addClass(getOut);
                        setTimeout(function() {
                            $(getId).removeClass("in");
                            $(getId).removeClass(getIn)
                        }, 500);
                        return false
                    })
                })
            }
        },
        navbarSticky: function() {
            var getNav = $("nav.navbar.bootsnav"),
                navSticky = getNav.hasClass("navbar-sticky");
            if(navSticky) {
                getNav.wrap("<div class='wrap-sticky'></div>");
                var getHeight = getNav.height();
                $(".wrap-sticky").height(getHeight);
                var getOffset = $(".wrap-sticky").offset().top;
                $(window).on("scroll", function() {
                    var scrollTop = $(window).scrollTop();
                    if(scrollTop > getOffset) {
                        getNav.addClass("sticked")
                    } else {
                        getNav.removeClass("sticked")
                    }
                })
            }
        },
        navbarScrollspy: function() {
            var navScrollSpy = $(".navbar-scrollspy"),
                $body = $('body'),
                getNav = $('nav.navbar.bootsnav'),
                offset = getNav.outerHeight();
            if(navScrollSpy.length) {
                $body.scrollspy({
                    target: '.navbar',
                    offset: offset
                });
                $('.scroll').on('click', function(event) {
                    event.preventDefault();
                    $('.scroll').removeClass("active");
                    $(this).addClass("active");
                    $(".navbar-collapse").removeClass("in");
                    $(".navbar-toggle").each(function() {
                        $(".fa", this).removeClass("fa-times");
                        $(".fa", this).addClass("fa-bars")
                    });
                    var scrollTop = $(window).scrollTop(),
                        $anchor = $(this).find('a'),
                        $section = $($anchor.attr('href')).offset().top,
                        $window = $(window).width(),
                        $minusDesktop = getNav.data("minus-value-desktop"),
                        $minusMobile = getNav.data("minus-value-mobile"),
                        $speed = getNav.data("speed");
                    if($window > 992) {
                        var $position = $section - $minusDesktop
                    } else {
                        var $position = $section - $minusMobile
                    }
                    $('html, body').stop().animate({
                        scrollTop: $position
                    }, $speed)
                });
                var fixSpy = function() {
                    var data = $body.data('bs.scrollspy');
                    if(data) {
                        offset = getNav.outerHeight();
                        data.options.offset = offset;
                        $body.data('bs.scrollspy', data);
                        $body.scrollspy('refresh')
                    }
                }
                var resizeTimer;
                $(window).on('resize', function() {
                    clearTimeout(resizeTimer);
                    var resizeTimer = setTimeout(fixSpy, 200)
                })
            }
        }
    };
    $(document).ready(function() {
        bootsnav.initialize()
    });
    $(window).on("resize", function() {
        bootsnav.hoverDropdown();
        bootsnav.navbarSticky();
        $(".navbar-toggle").each(function() {
            $(".fa", this).removeClass("fa-times");
            $(".fa", this).addClass("fa-bars");
            $(this).removeClass("fixed")
        });
        $(".navbar-collapse").removeClass("in");
        $(".navbar-collapse").removeClass("on");
        $(".navbar-collapse").removeClass("bounceIn")
    })
}(jQuery));