# -*- coding: utf-8 -*-

# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html

import pandas as pd


class HuaweirankPipeline(object):
    def open_spider(self, spider):
        self.df_result = []

    def process_item(self, item, spider):
        self.df_result.append(item)
        return item

    def close_spider(self, spider):
        result = pd.DataFrame(self.df_result)
        result.to_csv("rank.csv", index=False, encoding="utf_8_sig")
