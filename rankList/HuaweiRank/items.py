# -*- coding: utf-8 -*-

# Define here the models for your scraped items
#
# See documentation in:
# https://docs.scrapy.org/en/latest/topics/items.html

import scrapy


class HuaweirankItem(scrapy.Item):
    # define the fields for your item here like:
    # name = scrapy.Field()
    division = scrapy.Field()
    score = scrapy.Field()
    # rank = scrapy.Field()
    team_name = scrapy.Field()
    users = scrapy.Field()
    submit_time = scrapy.Field()
